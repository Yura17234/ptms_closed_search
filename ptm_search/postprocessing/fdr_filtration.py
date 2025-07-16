'''
    Фильтрация результатов PTM поиска методом target-decoy
'''

import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error

def threshold_calculation_identipy(df_Decoy, df_Target, log_file):
    FDR_list = []
    thresholds_q_values_dict = {}

    df_Decoy['decoy'] = [True] * df_Decoy.shape[0]
    df_Target['decoy'] = [False] * df_Target.shape[0]
    sum_df = pd.concat([df_Decoy, df_Target], ignore_index=True)
    sum_df = sum_df.sort_values(by=["log_hyperscore"], ascending=True)
    sum_df['rank'] = [k for k in range(1, len(sum_df) + 1)]
    df_Decoy = sum_df.query('decoy == True')
    df_Target = sum_df.query('decoy == False')

    for i in tqdm(np.linspace(int(df_Decoy['rank'].min()), int(df_Decoy['rank'].max()), 100000, dtype=int)[::-1]):
        if df_Target.query(f'rank >= {i}').shape[0] == 0:
            print('BAD')
            print(f'FDR: {FDR_list[-1]}, rank threshold: {i}')
            break
        FDR = df_Decoy.query(f'rank >= {i}').shape[0] / ( df_Target.query(f'rank >= {i}').shape[0] )
        FDR_list.append(FDR)
        FDR_threshold = i
        thresholds_q_values_dict[i] = FDR
        if FDR > 0.01:
            if round(FDR_list[-1], 2) <= 0.01:
                print(f'rounded FDR value: {round(FDR_list[-1], 2)}')
                log_file.write(f'rounded FDR value: {round(FDR_list[-1], 2)}\n')
                print('===============')
                print(FDR_list[-1], FDR_threshold)
                log_file.write(f'===============\nFDR: {FDR_list[-1]}, rank threshold: {FDR_threshold}\n\n')
                return FDR_threshold, thresholds_q_values_dict
            print('BAD')
            print(FDR_list[-1], FDR_threshold)
            log_file.write(f'BAD\n{FDR_list[-1]}, {FDR_threshold}\n\n')
            break

        if FDR <= 0.01 and FDR >= 0.0095:
            print('===============')
            print(f'FDR: {FDR}, rank threshold: {i}')
            log_file.write(f'\n===============\nFDR: {FDR}, rank threshold: {i}')
            return FDR_threshold, thresholds_q_values_dict

def calc_err_threshold(thresholds, errs):
    mean_error = np.mean(errs)
    std_error = np.std(errs)
    z_score_threshold = mean_error + 1 * std_error  # можно 3*std — для более строгого отсева

    # Находим первый индекс, где ошибка превышает порог
    idx1 = np.argmax(errs > z_score_threshold)
    err_threshold = thresholds[idx1] if idx1 < len(thresholds) else None
    return err_threshold

# ----------------------------------------------------------------------------------------------------------------------
def generate_knots(first_ref_threshold_, rank_before_err_, index_before_err_, thresholds_):
    if rank_before_err_ != 0:
        thresholds_before_err_ = thresholds_[0:index_before_err_]
        return thresholds_before_err_, np.sort(np.array(
            np.linspace(0, first_ref_threshold_, 10, dtype=int).tolist() +
            [int(round((first_ref_threshold_ + rank_before_err_) / 2, 0)), rank_before_err_]
        )).reshape(-1, 1)
    elif rank_before_err_ != 0 and rank_before_err_ < first_ref_threshold_:
        thresholds_before_err_ = thresholds_[0:index_before_err_]
        return thresholds_before_err_, np.sort(np.array(
        np.linspace(0, rank_before_err_, 12, dtype=int).tolist()
        )).reshape(-1, 1)
    else:
        thresholds_before_err_ = thresholds_[0:index_before_err_]
        return thresholds_before_err_, np.sort(np.array(
            np.linspace(0, first_ref_threshold_, 6, dtype=int).tolist() +
            np.linspace(first_ref_threshold_, max(thresholds_before_err_), 5, dtype=int).tolist()[1:]
        )).reshape(-1, 1)

def threshold_calculation_for_PTM_by_ranks(df_decoy_ss_and_ptm, df_target_ss_and_ptm, log_dir, log_file, config, ptm_name):
    # df_target_ss = df_target_ss_and_ptm.query("PTM == '-'")
    df_target_ptm = df_target_ss_and_ptm.query("PTM == '+'")
    # df_decoy_ss = df_decoy_ss_and_ptm.query("PTM == '-'")
    df_decoy_ptm = df_decoy_ss_and_ptm.query("PTM == '+'")

    # ------------------------------------------------------------------------------------------------------------------
    proportions_decoy_target, thresholds_decoy_target = [], []

    thresholds_decoy_target = np.linspace(0, int(df_target_ss_and_ptm['rank'].max()), 500, dtype=int)
    for q in tqdm(thresholds_decoy_target):
        x = df_decoy_ss_and_ptm.query(f'rank >= {q}').shape[0]
        y = df_target_ss_and_ptm.query(f'rank >= {q}').shape[0]

        if y == 0:
            break
        proportions_decoy_target.append(x / y)

    sigma_val = len(np.gradient(proportions_decoy_target)) / 100
    gaussian_filter_proportions_decoy_target = gaussian_filter(np.gradient(proportions_decoy_target), sigma=sigma_val)

    index_of_min = np.where(
        np.min(gaussian_filter_proportions_decoy_target) == gaussian_filter_proportions_decoy_target)
    first_ref_threshold = thresholds_decoy_target[index_of_min][0]
    # ------------------------------------------------------------------------------------------------------------------
    thresholds = np.linspace(df_decoy_ptm['rank'].min(), df_decoy_ptm['rank'].max(), 1000, dtype=int)

    proportions, error_propagation = [], []
    rank_before_err, index_before_err = 0, 0
    for index, q in tqdm(enumerate(thresholds)):
        x = df_decoy_ptm.query(f'rank >= {q}').shape[0]
        y = df_decoy_ss_and_ptm.query(f'rank >= {q}').shape[0]

        if y == 0:
            break

        err = (x / y) * np.sqrt((np.sqrt(x) / x) ** 2 + (np.sqrt(y) / y) ** 2)
        error_propagation.append(err)
        proportions.append(x / y)
        if err > 0.01 and rank_before_err == 0:
            rank_before_err = thresholds[index - 1]
            index_before_err = index - 1
        if rank_before_err == 0:
            index_before_err = index

    # ------------------------------------------------------------------------------------------------------------------
    errors = np.array(error_propagation)
    err_thresholds = np.array(thresholds)

    threshold1 = calc_err_threshold(err_thresholds, errors)

    err_sigma = len(np.gradient(errors)) / 10
    gaussian_error_propagation = gaussian_filter(np.gradient(errors), sigma=err_sigma)
    gaus_threshold = calc_err_threshold(err_thresholds, gaussian_error_propagation)

    # ------------------------------------------------------------------------------------------------------------------

    thresholds_before_err, knots = generate_knots(first_ref_threshold, rank_before_err, index_before_err, thresholds)

    X_train, X_test, y_train, y_test = train_test_split(np.array(thresholds_before_err).reshape((-1, 1)),
                                                        np.array(proportions[0:index_before_err]),
                                                        test_size=0.3, random_state=42)
    spl = make_pipeline(
        SplineTransformer(degree=1, knots=knots, extrapolation='linear'),
        Ridge()
    )
    spl.fit(X_train, y_train)
    y_pred = spl.predict(np.array(thresholds).reshape((-1, 1)))
    df_spline = pd.DataFrame({'thresholds': [float(w) for w in np.array(thresholds).reshape((-1, 1))],
                              'predicted_proportions': [float(w) for w in y_pred]}).sort_values(by=["thresholds"],
                                                                                                ascending=False)

    # ------------------------------------------------------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 15))
    sns.set(font_scale=1)
    sns.set_style("darkgrid", {"grid.color": ".1", "grid.linestyle": ":"})

    plt.subplot(2, 2, 1).set_title('Decoy / Target')
    plt.scatter(thresholds_decoy_target, proportions_decoy_target, color='#ff6666', s=40, alpha=0.4,
                edgecolors="#404040")
    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    plt.xlabel("Порог rank")
    plt.ylabel("Пропорция (Decoy / Target)")

    plt.subplot(2, 2, 2).set_title('Производная ( Decoy / Target )')
    plt.plot(thresholds_decoy_target, np.gradient(proportions_decoy_target))
    plt.plot(thresholds_decoy_target, gaussian_filter_proportions_decoy_target)
    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    plt.xlabel("Порог rank")
    plt.ylabel("Производная пропорции (Decoy / Target) ")

    plt.subplot(2, 2, 3).set_title('Пропорция decoy (PTM / Unmodified)')  # proportions approximation
    plt.scatter(thresholds, proportions, color='#3399ff', s=40, alpha=0.4, edgecolors="#404040", label='Пропорция')
    plt.plot(df_spline['thresholds'], df_spline['predicted_proportions'], color='red', label='Сплайн регрессии')

    r_squared = r'$R^2 =$' + str(
        round(r2_score(proportions[0:index_before_err], spl.predict(np.array(thresholds_before_err).reshape((-1, 1)))),
              2))
    # log_file.write(f'R^2 = {round(r2_score(proportions[0:index_before_err], spl.predict(np.array(thresholds_before_err).reshape((-1, 1)))), 2)}\n')
    plt.text(X_test.max() * (85 / 100), y_test.max() * (90 / 100), r_squared, weight='bold',
             horizontalalignment='center')
    rmse = r'$RMSE =$' + str(round(np.sqrt(mean_squared_error(proportions[0:index_before_err], spl.predict(
        np.array(thresholds_before_err).reshape((-1, 1))))), 4))
    # log_file.write(f'RMSE = {round(np.sqrt(mean_squared_error(proportions[0:index_before_err], spl.predict(np.array(thresholds[0:index_before_err]).reshape((-1, 1))))), 4)}\n')
    plt.text(X_test.max() * (85 / 100), y_test.max() * (85 / 100), rmse, weight='bold', horizontalalignment='center')

    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    if rank_before_err != 0:
        plt.axvline(rank_before_err, color="#404040", linestyle="-")
    plt.axvline(threshold1, color="#FF7B00", linestyle="-")
    plt.axvline(gaus_threshold, color="#0044FF", linestyle="-")

    plt.xlabel('Порог rank')  # threshhold rank
    plt.ylabel('Пропорция')  # proportion
    plt.legend()

    plt.subplot(2, 2, 4).set_title('Ошибка частного')  # error propagation
    plt.scatter(thresholds, error_propagation, color='#33cc33', s=40, alpha=0.4, edgecolors="#404040")
    if rank_before_err != 0:
        plt.axvline(rank_before_err, color="#404040", linestyle="-")
    plt.axvline(threshold1, color="#FF7B00", linestyle="-")
    plt.axvline(gaus_threshold, color="#0044FF", linestyle="-")

    plt.xlabel("Порог rank")  # threshold rank
    plt.ylabel("Ошибка частного")  # error propagation

    plt.figtext(0.5, 0.9, f'{ptm_name}', ha='center', va='center')

    plt.savefig(log_dir / f"{ptm_name.replace(' ', '_')}_Proportion_and_spline_regression.png", dpi=100,
                bbox_inches='tight')
    # ------------------------------------------------------------------------------------------------------------------
    # Вычисление попрога FDR на уровне 1% для PTM идентификаций
    fdr_threshold, fdrs_ptm_list, thresholds_q_values_dict = 0, [], {}
    for i in tqdm(thresholds[::-1]):
        fdr = df_decoy_ss_and_ptm.query(f'rank >= {i}').shape[0] / df_target_ss_and_ptm.query(f'rank >= {i}').shape[0]

        try:
            lambda_coef = df_target_ss_and_ptm.query(f'rank >= {i}').shape[0] / \
                          df_target_ptm.query(f'rank >= {i}').shape[0]
            gamma_coef = spl.predict(np.array([i]).reshape((-1, 1)))[0]
            fdr_ptm = lambda_coef * gamma_coef * fdr
            fdrs_ptm_list.append(fdr_ptm)
        except:
            print('BAD')
            # print(fdrs_ptm_list[-1], fdr_threshold)
            return fdr_threshold, thresholds_q_values_dict

        # print(FDR_PTM, FDR_threshold)
        fdr_threshold = i
        thresholds_q_values_dict[i] = fdr_ptm
        if fdr_ptm > 0.01:
            if round(fdrs_ptm_list[-1], 2) <= 0.01:
                print(f'rounded FDR value: {round(fdrs_ptm_list[-1], 2)}')
                # log_file.write(f'rounded FDR value: {round(fdrs_ptm_list[-1], 2)}\n')
                print('===============')
                print(fdrs_ptm_list[-1], fdr_threshold)
                # log_file.write(f'===============\nFDR: {fdrs_ptm_list[-1]}, rank threshold: {fdr_threshold}\n\n')
                return fdr_threshold, thresholds_q_values_dict
            print('BAD')
            print(fdrs_ptm_list[-1], fdr_threshold)
            return fdr_threshold, thresholds_q_values_dict

            # log_file.write(f'BAD\n{fdrs_ptm_list[-1]}, {fdr_threshold}\n\n')
            # break

        if fdr_ptm <= 0.01 and fdr_ptm >= 0.0095:  # 0.0089 | 0.005
            print('===============')
            print(fdr_ptm, fdr_threshold)
            # log_file.write(f'===============\nFDR: {fdr_ptm}, rank threshold: {fdr_threshold}\n\n')
            return fdr_threshold, thresholds_q_values_dict
