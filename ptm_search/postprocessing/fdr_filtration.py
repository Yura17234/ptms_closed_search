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
from sklearn.linear_model import RidgeCV
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

def threshold_calculation_for_PTM_by_ranks(df_Decoy_FS_and_PTM, df_Target_FS_and_PTM, log_dir, log_file, config, ptm_name):
    df_Target_FS = df_Target_FS_and_PTM.query("PTM == '-'")
    df_Target_PTM = df_Target_FS_and_PTM.query("PTM == '+'")
    df_Decoy_FS = df_Decoy_FS_and_PTM.query("PTM == '-'")
    df_Decoy_PTM = df_Decoy_FS_and_PTM.query("PTM == '+'")

    # ------------------------------------------------------------------------------------------------------------------
    proportions_decoy_target = []
    thresholds_decoy_target = []

    for q in tqdm(np.linspace(0, int(df_Target_FS['rank'].max()), 1000, dtype=int)):
        if df_Target_PTM.query(f'rank >= {q}').shape[0] == 0 and df_Target_FS.query(f'rank >= {q}').shape[0] == 0:
            break
        x = df_Decoy_PTM.query(f'rank >= {q}').shape[0] + df_Decoy_FS.query(f'rank >= {q}').shape[0]
        y = df_Target_PTM.query(f'rank >= {q}').shape[0] + df_Target_FS.query(f'rank >= {q}').shape[0]

        proportions_decoy_target.append(x / y)
        thresholds_decoy_target.append(q)

    gaussian_filter_proportions_decoy_target = gaussian_filter(np.gradient(proportions_decoy_target), sigma=20)
    first_threshold_index = [x[0] for x in zip(thresholds_decoy_target, gaussian_filter_proportions_decoy_target) if x[1] == min(gaussian_filter_proportions_decoy_target)][0]
    # for u in range(0, gaussian_filter_proportions_decoy_target.shape[0]):
    #     if gaussian_filter_proportions_decoy_target[u] > gaussian_filter_proportions_decoy_target[u - 1]:
    #
    #         first_threshold_index = thresholds_decoy_target[u]
    #         break
    # ------------------------------------------------------------------------------------------------------------------
    proportions = []
    thresholds = []
    error_propagation = []

    for q in tqdm(np.linspace(int(df_Decoy_PTM['rank'].min()), int(df_Decoy_PTM['rank'].max()), 1000, dtype=int)):
        if df_Decoy_PTM.query(f'rank >= {q}').shape[0] == 0 and df_Decoy_FS.query(f'rank >= {q}').shape[0] == 0:
            break

        x = df_Decoy_PTM.query(f'rank >= {q}').shape[0]
        y = df_Decoy_PTM.query(f'rank >= {q}').shape[0] + df_Decoy_FS.query(f'rank >= {q}').shape[0]

        error_propagation.append( (x / y) * np.sqrt((np.sqrt(x) / x) ** 2 + (np.sqrt(y) / y) ** 2) )
        proportions.append(x / y)
        thresholds.append(q)

    last_rank = 0
    index = 0
    for u in zip(thresholds, error_propagation):
        if u[1] >= 0.01 and last_rank == 0:
            last_rank = thresholds[index - 1]
            break
        index += 1
    # print(error_propagation)
    # print(last_rank, index)
    # ------------------------------------------------------------------------------------------------------------------
    knots = []
    if last_rank != 0:
        knots = np.linspace(0, first_threshold_index, 10, dtype=int).tolist() + [int(round((first_threshold_index+last_rank)/2, 0)), last_rank]
        knots = np.array(knots).reshape(-1, 1)
    elif last_rank == 0:
        knots = np.linspace(0, first_threshold_index, 6, dtype=int).tolist() + np.linspace(first_threshold_index, max(thresholds[0:index]), 5, dtype=int).tolist()[1:]
        knots = np.array(knots).reshape(-1, 1)

    X_train, X_test, y_train, y_test = train_test_split(np.array(thresholds[0:index]).reshape((-1, 1)), np.array(proportions[0:index]),
                                                        test_size=0.3, random_state=42)

    spl = make_pipeline(
        SplineTransformer(degree=1, knots=knots, extrapolation='linear'),
        RidgeCV(alphas=0.01)
    )
    spl.fit(X_train, y_train)
    y_pred = spl.predict(np.array(thresholds).reshape((-1, 1)))

    df_spline = pd.DataFrame({'thresholds': [float(w) for w in np.array(thresholds).reshape((-1, 1))],
                       'predicted_proportions': [float(w) for w in y_pred]}).sort_values(by=["thresholds"], ascending=False)

    # ------------------------------------------------------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 15))
    sns.set(font_scale=1)
    sns.set_style("darkgrid", {"grid.color": ".1", "grid.linestyle": ":"})

    plt.subplot(2, 2, 1).set_title('Decoy / Target')
    plt.scatter(thresholds_decoy_target, proportions_decoy_target, color='#ff6666', s=40, alpha=0.4,
                edgecolors="#404040")
    plt.axvline(x=first_threshold_index, color='red', linestyle="--")
    plt.xlabel("Порог rank")
    plt.ylabel("Пропорция (Decoy / Target)")

    plt.subplot(2, 2, 2).set_title('Производная ( Decoy / Target )')
    plt.plot(thresholds_decoy_target, np.gradient(proportions_decoy_target))
    plt.plot(thresholds_decoy_target, gaussian_filter_proportions_decoy_target)
    plt.axvline(x=first_threshold_index, color='red', linestyle="--")
    plt.xlabel("Порог rank")
    plt.ylabel("Производная пропорции (Decoy / Target) ")

    plt.subplot(2, 2, 3).set_title('Пропорция decoy (PTM / Unmodified)') # proportions approximation
    plt.scatter(thresholds, proportions, color='#3399ff', s=40, alpha=0.4, edgecolors="#404040", label='Пропорция')
    plt.plot(df_spline['thresholds'], df_spline['predicted_proportions'], color='red', label='Сплайн регрессии')

    r_squared = r'$R^2 =$' + str(round(r2_score(proportions[0:index], spl.predict(np.array(thresholds[0:index]).reshape((-1, 1)))), 2))
    log_file.write(f'R^2 = {round(r2_score(proportions[0:index], spl.predict(np.array(thresholds[0:index]).reshape((-1, 1)))), 2)}\n')
    plt.text(X_test.max()*(85/100), y_test.max()*(90/100), r_squared, weight='bold', horizontalalignment='center')
    rmse = r'$RMSE =$' + str(round(np.sqrt(mean_squared_error(proportions[0:index], spl.predict(np.array(thresholds[0:index]).reshape((-1, 1))))), 4))
    log_file.write(f'RMSE = {round(np.sqrt(mean_squared_error(proportions[0:index], spl.predict(np.array(thresholds[0:index]).reshape((-1, 1))))), 4)}\n')
    plt.text(X_test.max() * (85 / 100), y_test.max() * (85 / 100), rmse, weight='bold', horizontalalignment='center')

    plt.axvline(x=first_threshold_index, color='red', linestyle="--")
    if last_rank != 0:
        plt.axvline(last_rank, color="#404040", linestyle="-")
    plt.xlabel('Порог rank') # threshhold rank
    plt.ylabel('Пропорция') # proportion
    plt.legend()

    plt.subplot(2, 2, 4).set_title('Ошибка частного') # error propagation
    plt.scatter(thresholds, error_propagation, color='#33cc33', s=40, alpha=0.4, edgecolors="#404040")
    if last_rank != 0:
        plt.axvline(last_rank, color="#404040", linestyle="-")
    plt.xlabel("Порог rank") # threshold rank
    plt.ylabel("Ошибка частного") # error propagation

    plt.figtext(0.5, 0.9, f'{ptm_name}', ha='center', va='center')
    plt.savefig(log_dir / f"{ptm_name.replace(' ', '_')}_Proportion_and_spline_regression.png", dpi=100,
                bbox_inches='tight')
    plt.close(fig)
    # ------------------------------------------------------------------------------------------------------------------
    # Вычисление попрога FDR на уровне 1% для PTM идентификаций
    FDR_threshold = 0
    list_FDRs_PTM = []
    thresholds_q_values_dict = {}
    for i in tqdm(np.linspace(int(df_Decoy_PTM['rank'].min()), int(df_Decoy_PTM['rank'].max()), 1000, dtype=int)[::-1]): # [::-1]
        FDR_FS_and_PTM = df_Decoy_FS_and_PTM.query(f'rank >= {i}').shape[0] / ( df_Target_FS_and_PTM.query(f'rank >= {i}').shape[0] )# + df_Decoy_FS_and_PTM.query(f'rank >= {i}').shape[0] )

        try:
            FDR_PTM = (( df_Target_FS.query(f'rank >= {i}').shape[0]+df_Target_PTM.query(f'rank >= {i}').shape[0] ) / df_Target_PTM.query(f'rank >= {i}').shape[0]) * (spl.predict(np.array([i]).reshape((-1, 1)))[0]) * (FDR_FS_and_PTM)
            list_FDRs_PTM.append(FDR_PTM)
        except:
            print('BAD')
            try:
                print(list_FDRs_PTM[-1], FDR_threshold)
                log_file.write(f'BAD\n{list_FDRs_PTM[-1]}, {FDR_threshold}\n\n')
            except:
                break
            break

        # print(FDR_PTM, FDR_threshold)
        FDR_threshold = i
        thresholds_q_values_dict[i] = FDR_PTM
        if FDR_PTM > 0.01:
            if round(list_FDRs_PTM[-1], 2) <= 0.01:
                print(f'rounded FDR value: {round(list_FDRs_PTM[-1], 2)}')
                log_file.write(f'rounded FDR value: {round(list_FDRs_PTM[-1], 2)}\n')
                print('===============')
                print(list_FDRs_PTM[-1], FDR_threshold)
                log_file.write(f'===============\nFDR: {list_FDRs_PTM[-1]}, rank threshold: {FDR_threshold}\n\n')
                return FDR_threshold, thresholds_q_values_dict
            print('BAD')
            print(list_FDRs_PTM[-1], FDR_threshold)
            log_file.write(f'BAD\n{list_FDRs_PTM[-1]}, {FDR_threshold}\n\n')
            break

        if FDR_PTM <= 0.01 and FDR_PTM >= 0.0095:# 0.0089 | 0.005
            print('===============')
            print(FDR_PTM, FDR_threshold)
            log_file.write(f'===============\nFDR: {FDR_PTM}, rank threshold: {FDR_threshold}\n\n')
            return FDR_threshold, thresholds_q_values_dict
