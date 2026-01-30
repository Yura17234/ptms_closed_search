import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
from pyteomics.auxiliary.target_decoy import filter
from dataclasses import dataclass
from pathlib import Path
import logging

logger = logging.getLogger("aggregate_results")

'''
    Filtering PTM search results using the target-decoy approach
'''

def threshold_calculation_separate_fdr_by_pyteomics(df_decoy: pd.DataFrame, df_target: pd.DataFrame) -> pd.DataFrame:
    df_psms = pd.concat([df_decoy, df_target], ignore_index=True)
    df_psms['expect'] = 1 / df_psms['hyperscore']
    return filter(df_psms, key='expect', is_decoy=df_psms['decoy'], fdr=0.01, desc=True)

def threshold_calculation_separate_fdr(df_decoy: pd.DataFrame, df_target: pd.DataFrame, fdr_level=0.01) -> tuple[float, dict[float, float]]:

    decoys = np.sort(df_decoy['rank'].to_numpy())
    targets = np.sort(df_target['rank'].to_numpy())

    thresholds = np.sort(np.concatenate([df_decoy['rank'].to_numpy(), df_target['rank'].to_numpy()]))

    fdr_threshold, fdr_list, thresholds_q_values_dict = 0, [], {}
    for t in tqdm(thresholds):

        n_decoy = len(decoys) - np.searchsorted(decoys, t, side="left")
        n_target = len(targets) - np.searchsorted(targets, t, side="left")

        if n_target == 0:
            logger.info(f'BAD -- FDR: {fdr_list[-1]:.5f}, rank threshold: {t}')
            continue

        fdr = n_decoy / n_target

        fdr_threshold = t
        fdr_list.append(fdr)
        thresholds_q_values_dict[t] = fdr

        if fdr <= fdr_level:
            logger.info('---------------')
            logger.info(f'FDR : {fdr:.5f}, rank threshold: {t}')
            return t, thresholds_q_values_dict

    if fdr_list:
        logger.info(f"BAD -- min FDR PTM: {fdr_list[-1]:.5f}, Threshold: {fdr_threshold}")
    else:
        logger.info("BAD -- No valid FDR values found.")

@dataclass
class ApproximationPlotPayload:
    ratio_info_dir: Path
    ptm_name: str
    thresholds_decoy_target: np.ndarray
    proportions_decoy_target: np.ndarray
    gaussian_filter_proportions_decoy_target: np.ndarray
    thresholds: np.ndarray
    ptm_proportions: np.ndarray
    error_propagation: np.ndarray
    df_spline: pd.DataFrame
    first_ref_threshold: int
    index_before_err: int
    thresholds_before_err: np.ndarray
    rank_before_err: int

def plot_approximation(payload, spl) -> None:
    sns.set(font_scale=1)
    sns.set_style('darkgrid', {'grid.color': '.1', 'grid.linestyle': ':'})

    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    ax1, ax2, ax3, ax4 = axes.ravel()

    def vlines(ax):
        ax.axvline(x=payload.first_ref_threshold, color='red', linestyle='--')
        if getattr(payload, 'rank_before_err', 0):
            ax.axvline(x=payload.rank_before_err, color='green', linestyle='--')

    # --- (1) Proportion Decoy/Target ---
    ax1.set_title('Proportion (Decoy / Target)')
    ax1.scatter(payload.thresholds_decoy_target, payload.proportions_decoy_target,
        color='#ff6666', s=40, alpha=0.4, edgecolors='#404040')
    vlines(ax1)
    ax1.set_xlabel('rank')
    ax1.set_ylabel('Proportion')

    # --- (2) Gradient + gaussian filtered gradient ---
    grad = np.gradient(payload.proportions_decoy_target)
    ax2.set_title('Gradient of proportion (Decoy / Target)')
    ax2.plot(payload.thresholds_decoy_target, grad)
    ax2.plot(payload.thresholds_decoy_target, payload.gaussian_filter_proportions_decoy_target)
    vlines(ax2)
    ax2.set_xlabel('rank')
    ax2.set_ylabel("Gradient of proportion")

    # --- (3) PTM proportion + spline prediction ---
    ax3.set_title('Proportion of decoys (PTM / Unmodified + PTM)')
    ax3.scatter(payload.thresholds, payload.ptm_proportions,
        color='#3399ff', s=40, alpha=0.4, edgecolors='#404040', label='Proportion')

    ax3.plot(payload.df_spline['thresholds'], payload.df_spline['predicted_proportions'],
        color='red', label='Spline regression')

    vlines(ax3)
    ax3.set_xlabel("rank")
    ax3.set_ylabel("Proportion")
    ax3.legend()

    # metrics: R^2 and RMSE (train/fit diagnostics on available region)
    y_true = payload.ptm_proportions[: payload.index_before_err + 1]
    y_pred = spl.predict(np.asarray(payload.thresholds_before_err).reshape(-1, 1))

    r2 = r2_score(y_true, y_pred)
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    logger.info(f'R² = {r2:.4f}')
    logger.info(f'RMSE = {rmse:.4f}')

    # Place metrics relative to current axes limits (more stable than X_test.max())
    x0, x1 = ax3.get_xlim()
    y0, y1 = ax3.get_ylim()
    ax3.text(x0 + 0.85 * (x1 - x0), y0 + 0.90 * (y1 - y0), fr'$R^2 = {r2:.4f}$', weight='bold', ha='center')
    ax3.text(x0 + 0.85 * (x1 - x0), y0 + 0.85 * (y1 - y0), fr'$RMSE = {rmse:.4f}$', weight='bold', ha='center')

    # --- (4) Error propagation ---
    ax4.set_title('Error propagation')
    ax4.scatter(payload.thresholds, payload.error_propagation,
        color='#33cc33', s=40, alpha=0.4, edgecolors='#404040')

    if getattr(payload, 'rank_before_err', 0):
        ax4.axvline(x=payload.rank_before_err, color='green', linestyle='--')
    ax4.set_xlabel('rank')
    ax4.set_ylabel('Error propagation')

    fig.suptitle(str(payload.ptm_name), y=0.98)
    fig.tight_layout()

    fig.savefig(payload.ratio_info_dir / f"{payload.ptm_name.replace(' ', '_')}_Proportion_and_spline_regression.png",
                dpi=100, bbox_inches='tight')
    plt.close(fig)

def adaptive_spline_number_knots(n_points: int, *, max_points: int = 1000, max_knots: int = 12, min_knots: int = 3) -> int:
    k = int(round(n_points * max_knots / max_points))
    return int(np.clip(k, min_knots, max_knots))

def generate_knots(first_ref_threshold: int, rank_before_err: int, index_before_err: int, thresholds: np.ndarray, n_knots: int,) -> tuple[np.ndarray, np.ndarray]:
    thresholds_before_err = thresholds[: index_before_err + 1]

    def _knots_from_values(values: np.ndarray) -> np.ndarray:
        vals = np.unique(values.astype(int))
        return vals.reshape(-1, 1)

    max_thr = int(thresholds_before_err.max())

    # 'rank_before_err' after the reference threshold
    if rank_before_err and rank_before_err > first_ref_threshold:
        base = np.linspace(0, first_ref_threshold, max(n_knots - 2, 1), dtype=int)
        mid = int(round((first_ref_threshold + rank_before_err) / 2))
        values = np.concatenate([base, np.array([mid, rank_before_err], dtype=int)])
        return thresholds_before_err, _knots_from_values(values)

    # rank_before_err before reference threshold
    if rank_before_err and rank_before_err < first_ref_threshold:
        values = np.linspace(0, rank_before_err, max(n_knots, 1), dtype=int)
        return thresholds_before_err, _knots_from_values(values)

    # no rank_before_err
    k1 = max(int(round(n_knots / 2)), 1)

    left = np.linspace(0, first_ref_threshold, k1, dtype=int)
    right = np.linspace(first_ref_threshold, max_thr, n_knots - k1 + 1, dtype=int)[1:]
    values = np.concatenate([left, right])
    return thresholds_before_err, _knots_from_values(values)


def threshold_calculation_transferred_fdr(df_decoy_ss_and_ptm: pd.DataFrame, df_target_ss_and_ptm: pd.DataFrame, ratio_info_dir: Path, ptm_name: str, fdr_level=0.01) -> tuple[float, dict[float, float]]:
    all_decoys = np.sort(df_decoy_ss_and_ptm['rank'].to_numpy())
    all_targets = np.sort(df_target_ss_and_ptm['rank'].to_numpy())

    ptm_decoys = np.sort(df_decoy_ss_and_ptm[df_decoy_ss_and_ptm['PTM'] == True]['rank'].to_numpy())
    ptm_targets = np.sort(df_target_ss_and_ptm[df_target_ss_and_ptm['PTM'] == True]['rank'].to_numpy())

    proportions_decoy_target = []
    thresholds_decoy_target = np.linspace(1, df_target_ss_and_ptm['rank'].max(), 500, dtype=int)
    for all_r in tqdm(thresholds_decoy_target):
        n_all_decoy = len(all_decoys) - np.searchsorted(all_decoys, all_r, side='left')
        n_all_target = len(all_targets) - np.searchsorted(all_targets, all_r, side='left')

        if n_all_target == 0:
            break
        proportions_decoy_target.append(n_all_decoy / n_all_target)

    sigma_val = len(proportions_decoy_target) / 100
    gaussian_filter_proportions_decoy_target = gaussian_filter(np.gradient(proportions_decoy_target), sigma=sigma_val)

    idx_min = int(np.argmin(gaussian_filter_proportions_decoy_target))
    first_ref_threshold = int(thresholds_decoy_target[idx_min])

    # ------------------------------------------------------------------------------------------------------------------
    ptm_proportions, error_propagation = [], []
    rank_before_err, index_before_err = 0, 0
    thresholds = np.linspace(all_decoys.min(), all_decoys.max(), 1000, dtype=int)
    for index, ptm_p in enumerate(thresholds):
        x = len(ptm_decoys) - np.searchsorted(ptm_decoys, ptm_p, side='left')
        y = len(all_decoys) - np.searchsorted(all_decoys, ptm_p, side='left')

        if y == 0:
            break

        err = (x / y) * np.sqrt(1 / x + 1 / y)

        error_propagation.append(err)
        ptm_proportions.append(x / y)
        if err >= 0.01 and rank_before_err == 0:
            rank_before_err = thresholds[index - 1]
            index_before_err = index - 1
        elif rank_before_err == 0:
            index_before_err = index

    if max(error_propagation) <= 0.01 and rank_before_err == 0:
        rank_before_err = thresholds[index_before_err]

    ptm_proportions = np.array(ptm_proportions)
    error_propagation = np.array(error_propagation)

    # ------------------------------------------------------------------------------------------------------------------
    number_knots = adaptive_spline_number_knots(len(thresholds[0:index_before_err + 1]))

    thresholds_before_err, knots = generate_knots(first_ref_threshold, rank_before_err, index_before_err, thresholds,
                                                  number_knots)

    X_train, X_test, y_train, y_test = train_test_split(thresholds_before_err.reshape((-1, 1)),
                                                        ptm_proportions[0:index_before_err + 1],
                                                        test_size=0.2, random_state=42)
    spl = make_pipeline(
        SplineTransformer(degree=1, knots=knots, extrapolation='linear'),
        LinearRegression())

    spl.fit(X_train, y_train)
    y_pred = spl.predict(thresholds.reshape((-1, 1))).ravel()
    df_spline = pd.DataFrame({'thresholds': thresholds.astype(float),
                              'predicted_proportions': y_pred.astype(float)}).sort_values(by=["thresholds"], ascending=False)

    # ------------------------------------------------------------------------------------------------------------------
    payload = ApproximationPlotPayload(
        ratio_info_dir=ratio_info_dir,
        ptm_name=ptm_name,
        thresholds_decoy_target=thresholds_decoy_target,
        proportions_decoy_target=np.array(proportions_decoy_target),
        gaussian_filter_proportions_decoy_target=gaussian_filter_proportions_decoy_target,
        thresholds=thresholds,
        ptm_proportions=ptm_proportions,
        error_propagation=error_propagation,
        df_spline=df_spline,
        first_ref_threshold=int(first_ref_threshold),
        index_before_err=int(index_before_err),
        thresholds_before_err=thresholds_before_err,
        rank_before_err=int(rank_before_err))

    plot_approximation(payload, spl)

    # ------------------------------------------------------------------------------------------------------------------
    # Calculation of the FDR threshold at 1% for PTM identifications
    fdr_threshold, fdrs_ptm_list, thresholds_q_values_dict = 0, [], {}
    for i in tqdm(thresholds):
        n_all_decoy = len(all_decoys) - np.searchsorted(all_decoys, i, side='left')
        n_all_target = len(all_targets) - np.searchsorted(all_targets, i, side='left')

        if n_all_target == 0:
            break

        fdr = n_all_decoy / n_all_target

        try:
            n_ptm_target = len(ptm_targets) - np.searchsorted(ptm_targets, i, side='left')
            if n_ptm_target == 0:
                break
            lambda_coef = n_all_target / n_ptm_target
            gamma_coef = spl.predict(np.array([i]).reshape((-1, 1)))[0]
            gamma_coef = max(0.0, min(1.0, gamma_coef))

            fdr_ptm = lambda_coef * gamma_coef * fdr

        except Exception:
            if fdrs_ptm_list:
                logger.info(f'BAD -- min FDR PTM: {fdrs_ptm_list[-1]:.5f}, Threshold: {fdr_threshold}')
            else:
                logger.info("BAD -- first FDR PTM failed")
            break

        fdr_threshold = i
        thresholds_q_values_dict[i] = fdr_ptm
        fdrs_ptm_list.append(fdr_ptm)

        if fdr_ptm <= fdr_level:
            logger.info('---------------')
            logger.info(f'FDR : {fdr_ptm:.5f}, rank threshold: {i}')
            return i, thresholds_q_values_dict

    if fdrs_ptm_list:
        logger.info(f"BAD -- min FDR PTM: {fdrs_ptm_list[-1]:.5f}, Threshold: {fdr_threshold}")
    else:
        logger.info("BAD -- No valid FDR values found.")


def threshold_calculation_transferred_fdr_linear_reg(df_decoy_ss_and_ptm: pd.DataFrame, df_target_ss_and_ptm: pd.DataFrame, ratio_info_dir: Path, ptm_name: str) -> tuple[float, dict[float, float]]:
    df_target_ptm = df_target_ss_and_ptm.query("PTM == '+'")
    df_decoy_ptm = df_decoy_ss_and_ptm.query("PTM == '+'")

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
    thresholds_before_err = thresholds[0:index_before_err + 1]
    X_train, X_test, y_train, y_test = train_test_split(np.array(thresholds_before_err).reshape((-1, 1)),
                                                        np.array(ptm_proportions[0:index_before_err + 1]),
                                                        test_size=0.2, random_state=42)

    linear_reg_model = LinearRegression()
    linear_reg_model.fit(X_train, y_train)
    linear_f = np.poly1d(np.array(list(linear_reg_model.coef_)[::-1] + [linear_reg_model.intercept_]))
    linear_y_pred = linear_f(thresholds)

    df_spline = pd.DataFrame({'thresholds': [float(w) for w in np.array(thresholds).reshape((-1, 1))],
                              'predicted_proportions': [float(w) for w in linear_y_pred]}).sort_values(by=["thresholds"],
                                                                                                ascending=False)
    # ------------------------------------------------------------------------------------------------------------------

    fig = plt.figure(figsize=(15, 15))
    sns.set(font_scale=1)
    sns.set_style("darkgrid", {"grid.color": ".1", "grid.linestyle": ":"})

    plt.subplot(1, 2, 1).set_title('Proportion of decoys (PTM / Unmodified + PTM)')
    plt.scatter(thresholds, ptm_proportions, color='#3399ff', s=40, alpha=0.4, edgecolors="#404040", label='Пропорция')
    plt.plot(df_spline['thresholds'], df_spline['predicted_proportions'], color='red', label='Сплайн регрессии')

    # ------------------------------------------------------------------------------------------------------------------
    y_true = ptm_proportions[:index_before_err + 1]
    y_pred = linear_f(thresholds_before_err)

    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))

    plt.text(X_test.max() * 0.85, y_test.max() * 0.90,
             fr'$R^2 = {r2:.4f}$', weight='bold', ha='center')
    plt.text(X_test.max() * 0.85, y_test.max() * 0.85,
             fr'$RMSE = {rmse:.4f}$', weight='bold', ha='center')

    logger.info(f'R² = {r2:.4f}')
    logger.info(f'RMSE = {rmse:.4f}')
    # ------------------------------------------------------------------------------------------------------------------

    if rank_before_err != 0:
        plt.axvline(rank_before_err, color="green", linestyle="--")
    plt.xlabel('Threshold rank')
    plt.ylabel('Proportion')
    plt.legend()

    plt.subplot(1, 2, 2).set_title('Error propagation')
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
    # Calculation of the FDR threshold at 1% for PTM identifications
    fdr_threshold, fdrs_ptm_list, thresholds_q_values_dict = 0, [], {}
    for i in tqdm(np.linspace(df_decoy_ptm['rank'].min(), df_decoy_ptm['rank'].max(), 1000, dtype=int)[::-1]):
        fdr = df_decoy_ss_and_ptm.query(f'rank >= {i}').shape[0] / df_target_ss_and_ptm.query(f'rank >= {i}').shape[0]

        try:
            lambda_coef = df_target_ss_and_ptm.query(f'rank >= {i}').shape[0] / df_target_ptm.query(f'rank >= {i}').shape[0]
            gamma_coef = linear_f(i)
            fdr_ptm = lambda_coef * gamma_coef * fdr
        except:
            logger.info(f'BAD -- min FDR PTM: {fdrs_ptm_list[-1]:.4f}, Threshold: {fdr_threshold}')
            break

        fdr_threshold = i
        thresholds_q_values_dict[i] = fdr_ptm
        fdrs_ptm_list.append(fdr_ptm)

        if fdr_ptm > 0.01:
            if round(fdrs_ptm_list[-1], 2) <= 0.01:
                logger.info(f'rounded FDR value: {fdrs_ptm_list[-1]:.4f}')
                logger.info('---------------')
                logger.info(f'FDR PTM: {fdrs_ptm_list[-1]:.4f}, Threshold: {fdr_threshold}')
                return fdr_threshold, thresholds_q_values_dict
            logger.info(f'BAD -- min FDR PTM: {fdrs_ptm_list[-1]:.4f}, Threshold: {fdr_threshold}')
            break

        if fdr_ptm <= 0.01 and fdr_ptm >= 0.0095:  # 0.0089 | 0.005
            logger.info('---------------')
            logger.info(f'FDR PTM: {fdr_ptm:.4f}, Threshold: {fdr_threshold}')
            return fdr_threshold, thresholds_q_values_dict
