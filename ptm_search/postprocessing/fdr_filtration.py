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
from sklearn.metrics import r2_score, mean_squared_error
from pathlib import Path
import logging
logger = logging.getLogger("aggregate_results")

def threshold_calculation_separate_fdr(df_decoy: pd.DataFrame, df_target: pd.DataFrame) -> tuple[float, dict[float, float]]:
    fdr_threshold, fdr_list, thresholds_q_values_dict = 0, [], {}

    for i in tqdm(np.linspace(int(df_decoy['rank'].min()), int(df_target['rank'].max()), 1000, dtype=int)[::-1]):
        x = df_decoy.query(f'rank >= {i}').shape[0]
        y = df_target.query(f'rank >= {i}').shape[0]
        if y == 0:
            logger.info('BAD')
            logger.info(f'FDR: {fdr_list[-1]}, rank threshold: {i}')
            break
        fdr = x / y

        # if fdr > 0.01:
        #     try:
        #         t = (0.01 - fdr_list[-1]) / (fdr - fdr_list[-1])
        #         fdr_threshold = fdr_threshold + t * (i - fdr_threshold)
        #
        #         thresholds_q_values_dict[fdr_threshold] = fdr_list[-1]
        #         # print('===============')
        #         # print(fdrs_ptm_list[-1], fdr_threshold)
        #         print(f'===============\nFDR: {fdr_list[-1]}, rank threshold: {fdr_threshold}\n\n')
        #         return fdr_threshold, thresholds_q_values_dict
        #     except:
        #         if round(fdr, 2) <= 0.01:
        #             # print(f'rounded FDR value: {round(fdr_ptm, 2)}')
        #             # print('===============')
        #             print(fdr, i)
        #             thresholds_q_values_dict[i] = fdr
        #             print(f'===============\nFDR: {fdr}, rank threshold: {i}\n\n')
        #             return i, thresholds_q_values_dict
        #         print('BAD')
        #         print(fdr, fdr_threshold)
        #         return

        fdr_threshold = i
        thresholds_q_values_dict[i] = fdr
        fdr_list.append(fdr)

        if fdr > 0.01:
            if round(fdr_list[-1], 2) <= 0.01:
                logger.info(f'rounded FDR value: {round(fdr_list[-1], 2)}')
                # log_file.write(f'rounded FDR value: {round(fdr_list[-1], 2)}\n')
                logger.info('===============')
                logger.info(fdr_list[-1], fdr_threshold)
                # log_file.write(f'===============\nFDR: {fdr_list[-1]}, rank threshold: {fdr_threshold}\n\n')
                return fdr_threshold, thresholds_q_values_dict
            logger.info('BAD')
            logger.info(fdr_list[-1], fdr_threshold)
            # log_file.write(f'BAD\n{fdr_list[-1]}, {fdr_threshold}\n\n')
            break

        if fdr <= 0.01 and fdr >= 0.0095:
            logger.info('===============')
            logger.info(f'FDR: {fdr}, rank threshold: {i}')
            # log_file.write(f'\n===============\nFDR: {fdr}, rank threshold: {i}')
            return fdr_threshold, thresholds_q_values_dict

def adaptive_spline_number_knots(n_points, max_points=1000, max_knots=12, min_knots=3):
    knots = round(n_points * max_knots / max_points)
    return max(min_knots, min(knots, max_knots))

def generate_knots(first_ref_threshold_, rank_before_err_, index_before_err_, thresholds_, number_knots_):

    number_knots_1 = round(number_knots_ / 2)
    number_knots_2 = number_knots_1
    thresholds_before_err_ = thresholds_[0:index_before_err_+1]

    if rank_before_err_ != 0 and rank_before_err_ > first_ref_threshold_:
        return thresholds_before_err_, np.sort(np.array(
            np.linspace(0, first_ref_threshold_, number_knots_-2, dtype=int).tolist() +
            [int(round((first_ref_threshold_ + rank_before_err_) / 2, 0)), rank_before_err_]
        )).reshape(-1, 1)
    elif rank_before_err_ != 0 and rank_before_err_ < first_ref_threshold_:
        return thresholds_before_err_, np.sort(np.array(
        np.linspace(0, rank_before_err_, number_knots_, dtype=int).tolist()
        )).reshape(-1, 1)
    else:
        return thresholds_before_err_, np.sort(np.array(
            np.linspace(0, first_ref_threshold_, number_knots_1, dtype=int).tolist() +
            np.linspace(first_ref_threshold_, max(thresholds_before_err_), number_knots_2, dtype=int).tolist()[1:]
        )).reshape(-1, 1)

def threshold_calculation_transferred_fdr(df_decoy_ss_and_ptm: pd.DataFrame, df_target_ss_and_ptm: pd.DataFrame, ratio_info_dir: Path, ptm_name: str) -> tuple[float, dict[float, float]]:
    df_target_ptm = df_target_ss_and_ptm.query("PTM == '+'")
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

    ptm_proportions, error_propagation = [], []
    rank_before_err, index_before_err = 0, 0
    for index, q in tqdm(enumerate(thresholds)):
        x = df_decoy_ptm.query(f'rank >= {q}').shape[0]
        y = df_decoy_ss_and_ptm.query(f'rank >= {q}').shape[0]

        if y == 0:
            break

        err = (x / y) * np.sqrt((np.sqrt(x) / x) ** 2 + (np.sqrt(y) / y) ** 2)
        error_propagation.append(err)
        ptm_proportions.append(x / y)
        if err > 0.01 and rank_before_err == 0:
            rank_before_err = thresholds[index - 1]
            index_before_err = index - 1
        if rank_before_err == 0:
            index_before_err = index

    if max(error_propagation) <= 0.01 and rank_before_err == 0:
        rank_before_err = thresholds[index_before_err]

    # ------------------------------------------------------------------------------------------------------------------

    number_knots = adaptive_spline_number_knots(len(thresholds[0:index_before_err + 1]))

    thresholds_before_err, knots = generate_knots(first_ref_threshold, rank_before_err, index_before_err, thresholds,
                                                  number_knots)

    X_train, X_test, y_train, y_test = train_test_split(np.array(thresholds_before_err).reshape((-1, 1)),
                                                        np.array(ptm_proportions[0:index_before_err + 1]),
                                                        test_size=0.2, random_state=42)
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

    plt.subplot(2, 2, 1).set_title('Proportion (Decoy / Target)')
    plt.scatter(thresholds_decoy_target, proportions_decoy_target, color='#ff6666', s=40, alpha=0.4,
                edgecolors="#404040")
    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    plt.xlabel("Threshold rank")
    plt.ylabel("Proportion")

    plt.subplot(2, 2, 2).set_title('Gradient of proportion ( Decoy / Target )')
    plt.plot(thresholds_decoy_target, np.gradient(proportions_decoy_target))
    plt.plot(thresholds_decoy_target, gaussian_filter_proportions_decoy_target)
    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    plt.xlabel("Threshold rank")
    plt.ylabel("Gradient of proportion")

    plt.subplot(2, 2, 3).set_title('Proportion of decoys (PTM / Unmodified + PTM)')
    plt.scatter(thresholds, ptm_proportions, color='#3399ff', s=40, alpha=0.4, edgecolors="#404040", label='Пропорция')
    plt.plot(df_spline['thresholds'], df_spline['predicted_proportions'], color='red', label='Сплайн регрессии')

    # ------------------------------------------------------------------------------------------------------------------
    y_true = ptm_proportions[:index_before_err + 1]
    y_pred = spl.predict(np.array(thresholds_before_err).reshape(-1, 1))

    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))

    plt.text(X_test.max() * 0.85, y_test.max() * 0.90,
             fr'$R^2 = {r2:.4f}$', weight='bold', ha='center')
    plt.text(X_test.max() * 0.85, y_test.max() * 0.85,
             fr'$RMSE = {rmse:.4f}$', weight='bold', ha='center')

    logger.info(f'R² = {r2:.4f}')
    logger.info(f'RMSE = {rmse:.4f}\n')
    # ------------------------------------------------------------------------------------------------------------------

    plt.axvline(x=first_ref_threshold, color='red', linestyle="--")
    if rank_before_err != 0:
        plt.axvline(rank_before_err, color="green", linestyle="--")
    plt.xlabel('Threshold rank')
    plt.ylabel('Proportion')
    plt.legend()

    plt.subplot(2, 2, 4).set_title('Error propagation')
    plt.scatter(thresholds, error_propagation, color='#33cc33', s=40, alpha=0.4, edgecolors="#404040")
    if rank_before_err != 0:
        plt.axvline(rank_before_err, color="green", linestyle="--")
    plt.xlabel("Threshold rank")
    plt.ylabel("Error propagation")

    plt.figtext(0.5, 0.9, f'{ptm_name}', ha='center', va='center')
    plt.savefig(ratio_info_dir / f"{ptm_name.replace(' ', '_')}_Proportion_and_spline_regression.png", dpi=100,
                bbox_inches='tight')
    plt.close(fig)
    # ------------------------------------------------------------------------------------------------------------------
    # Вычисление попрога FDR на уровне 1% для PTM идентификаций
    fdr_threshold, fdrs_ptm_list, thresholds_q_values_dict = 0, [], {}
    for i in tqdm(np.linspace(df_decoy_ptm['rank'].min(), df_decoy_ptm['rank'].max(), 1000, dtype=int)[::-1]):
        fdr = df_decoy_ss_and_ptm.query(f'rank >= {i}').shape[0] / df_target_ss_and_ptm.query(f'rank >= {i}').shape[0]

        try:
            lambda_coef = df_target_ss_and_ptm.query(f'rank >= {i}').shape[0] / df_target_ptm.query(f'rank >= {i}').shape[0]
            gamma_coef = spl.predict(np.array([i]).reshape((-1, 1)))[0]
            fdr_ptm = lambda_coef * gamma_coef * fdr
        except:
            logger.info('BAD')
            logger.info(fdrs_ptm_list[-1], fdr_threshold)
            # return
            break

        # if fdr_ptm > 0.01:
        #     try:
        #         t = (0.01 - fdrs_ptm_list[-1]) / (fdr_ptm - fdrs_ptm_list[-1])
        #         fdr_threshold = fdr_threshold + t * (i - fdr_threshold)
        #
        #         thresholds_q_values_dict[fdr_threshold] = fdrs_ptm_list[-1]
        #         # print('===============')
        #         # print(fdrs_ptm_list[-1], fdr_threshold)
        #         print(f'===============\nFDR: {fdrs_ptm_list[-1]}, rank threshold: {fdr_threshold}\n\n')
        #         return fdr_threshold, thresholds_q_values_dict
        #     except:
        #         if round(fdr_ptm, 2) <= 0.01:
        #             # print(f'rounded FDR value: {round(fdr_ptm, 2)}')
        #             # print('===============')
        #             print(fdr_ptm, i)
        #             thresholds_q_values_dict[i] = fdr_ptm
        #             print(f'===============\nFDR: {fdr_ptm}, rank threshold: {i}\n\n')
        #             return i, thresholds_q_values_dict
        #         print('BAD')
        #         print(fdr_ptm, fdr_threshold)
        #         return

        fdr_threshold = i
        thresholds_q_values_dict[i] = fdr_ptm
        fdrs_ptm_list.append(fdr_ptm)


        if fdr_ptm > 0.01:
            if round(fdrs_ptm_list[-1], 2) <= 0.01:
                logger.info(f'rounded FDR value: {round(fdrs_ptm_list[-1], 2)}')
                # log_file.write(f'rounded FDR value: {round(fdrs_ptm_list[-1], 2)}\n')
                logger.info('===============')
                logger.info(fdrs_ptm_list[-1], fdr_threshold)
                # log_file.write(f'===============\nFDR: {fdrs_ptm_list[-1]}, rank threshold: {fdr_threshold}\n\n')
                return fdr_threshold, thresholds_q_values_dict
            logger.info('BAD')
            logger.info(fdrs_ptm_list[-1], fdr_threshold)
            # return fdr_threshold, thresholds_q_values_dict
            # return
            break

            # log_file.write(f'BAD\n{fdrs_ptm_list[-1]}, {fdr_threshold}\n\n')
            # break

        if fdr_ptm <= 0.01 and fdr_ptm >= 0.0095:  # 0.0089 | 0.005
            logger.info('===============')
            logger.info(fdr_ptm, fdr_threshold)
            # log_file.write(f'===============\nFDR: {fdr_ptm}, rank threshold: {fdr_threshold}\n\n')
            return fdr_threshold, thresholds_q_values_dict
