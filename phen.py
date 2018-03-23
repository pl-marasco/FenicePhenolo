#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import nodata
import outlier
import detect_peacks as dp
from seasonal import fit_seasons, fit_trend



def analyzes(dataset, prmts, coord):
    # Numpy error report
    np.seterr(all='ignore')  # {‘ignore’, ‘warn’, ‘raise’, ‘call’, ‘print’, ‘log’}

    # Error table

    err_table = pd.Series(0, index=['read', 'fix', 'rescaling', 'outlier', 'cleaning',
                                    's&t_ext', 's2d', 'savgol', 'valley'])

    # read x y
    try:
        h0 = dataset.read(prmts['bdsIdx'], window=coord)
    except RuntimeError:
        print('Error in read xy')
        err_table['read'] = 1
        return err_table

    coord = [coord[0][0], coord[1][0]]

    #
    # Core code equal across all the versions
    #

    # Water excluder
    if 254 in h0:
        return err_table
    else:
        h0 = h0.flatten().astype(np.float)

    # Output container
    ts_table = pd.DataFrame([])

    # Convert to pandas series
    ts_dek = pd.Series(h0, index=prmts['dates'])

    # Fix no data (works only for Spot Veg)
    try:
        ts_dek_unfix = ts_dek.copy()
        ts_dek = nodata.nodatafix(ts_dek)
    except (RuntimeError, Exception, ValueError):
        print('Error in fix no data, in position:{0}'.format(coord))
        err_table['fix'] = 1
        return err_table

    # Rescale values to  0-100
    try:
        # ts_dek_resc_unfix = (ts_dek_unfix - prmts['min_v']) / (prmts['max_v'] - prmts['min_v']) * 100
        ts_dek_resc = (ts_dek - prmts['min_v']) / (prmts['max_v'] - prmts['min_v']) * 100
    except (RuntimeError, Exception, ValueError):
        print('Error in rescaling, in position:{0}'.format(coord))
        err_table['rescaling'] = 1
        return err_table

    # Filter outlier
    try:
        ts_clean = outlier.season(ts_dek_resc, 1, prmts['dys_x_yr'] * prmts['outMax'], 3)  # TODO make variable dek
    except (RuntimeError, Exception, ValueError):
        print('Outlier research failed, in position:{0}'.format(coord))
        err_table['outlier'] = 1
        return err_table

    try:
        if ts_clean is not None:
            # ts_cleaned = outlier.fillerSeason(ts_clean)
            ts_cleaned = ts_clean.interpolate()
        else:
            # No season strategy
            summedYrs = ts_dek_resc.groupby(ts_dek_resc.index.year).sum()
            season_lng = 9999

            ts_table = {'sd': pd.Series(0, index=np.unique(ts_dek.index.year)),
                        'ed': pd.Series(0, index=np.unique(ts_dek.index.year)),
                        'sl': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
                        'spi': pd.Series(np.nan, index=np.unique(ts_dek.index.year)),
                        'si': summedYrs,
                        'cf': pd.Series(0, index=np.unique(ts_dek.index.year)),
                        'yr': pd.Series(np.unique(ts_dek.index.year), index=np.unique(ts_dek.index.year))}
            return pd.DataFrame(ts_table), err_table, season_lng
    except (RuntimeError, Exception, ValueError):
        print('Cleaning failed, in position:{0}'.format(coord))
        err_table['cleaning'] = 1
        return err_table

    # Estimate Season length
    try:
        seasons, trend = fit_seasons(ts_cleaned)
    except (RuntimeError, Exception, ValueError):
        print("Error in Season and Trend estimation, in position:{0}".format(coord))
        err_table['s&t_ext'] = 1
        return err_table

    # Calculate season length and expected number of season
    try:
        season_lng = len(seasons) * prmts['dys_mlt']
        expSeason = int((ts_cleaned.index.max() - ts_cleaned.index.min()) / pd.Timedelta(season_lng, unit='d'))
    except (RuntimeError, Exception, ValueError):
        print("Error! Season conversion to days failed, in position:{0}".format(coord))
        err_table['s2d'] = 1
        return err_table

    if prmts['medspan'] == 0:
        medspan = season_lng / 7
    else:
        medspan = prmts['medspan']

    # Interpolate data to daily sample
    ts_d = ts_cleaned.resample('D').asfreq().interpolate(method='linear').fillna(0)

    try:
        if prmts['smp'] != 0:
            # Savinsky Golet filter
            ps = savgol_filter(ts_d, prmts['medspan'], prmts['smp'], mode='nearest')
            # TODO automatic selection of savgol window
            ps = pd.Series(ps, ts_d.index)
        else:
            ps = ts_d.rolling(medspan // 2 * 2, win_type='boxcar', center=True).mean()
    except (RuntimeError, Exception, ValueError):
        print('Error! Savinsky Golet filter problem, in position:{0}'.format(coord))
        err_table['savgol'] = 1
        return err_table

    # Valley detection
    # Detrending to catch better points
    vtrend = pd.Series(fit_trend(ps), index=ps.index)
    vdetr = ps - vtrend
    try:
        if 200.0 < season_lng < 400.0:
            mpd_val = int(season_lng * 2/3)
        elif season_lng < 200:
            mpd_val = int(season_lng * 1/3)
        else:
            mpd_val = int(season_lng * (prmts['tr']-prmts['tr']*1/3) / 100)

        ind = dp.detect_peaks(vdetr, mph=vdetr.mean(),
                              mpd=mpd_val,
                              valley=True,

                              edge='both',
                              kpsh=False)
        oversmp = False

        if not ind.any():
            ind = dp.detect_peaks(vdetr, mph=-20,
                                  mpd=60,
                                  valley=True)
            oversmp = True

    except ValueError as e:
        print('Error in valley detection, in position:{0}, error {1}'.format(coord, e))
        err_table['vally'] = 1
        return err_table

    # Valley point time series conversion
    pks = pd.Series()
    for i in ind:
        pks[ps.index[i]] = ps.iloc[i]

    # region Plotting
    # plt.subplot(211)
    # mng = plt.get_current_fig_manager()
    # mng.window.state('zoomed')
    # ts_dek_unfix.plot(style='g--')
    # ts_dek.plot(style='r--')
    #
    # plt.subplot(212)
    # mng = plt.get_current_fig_manager()
    # mng.window.state('zoomed')
    # ps.plot(style='r:')
    # pks.plot(style='ro')
    # plt.show()
    # plt.clf()
    # pd.Series(-vdetr.values, index=vdetr.index).plot()
    # pd.Series((vdetr.std()), index=vdetr.index).plot()
    # end region

    pks0 = pd.Series()
    for i in ind:
        pks0[vdetr.index[i]] = -vdetr.iloc[i]

    pks0.plot(style='ro')
    # endregion

    # Point to point cycles
    for i in range(len(pks) - 1):

        # Min min curve
        mmi = ps[pks.index[i]: pks.index[i + 1]]

        # Interpolate line between min min
        pf = mmi.copy()
        pf.iloc[1:-1] = np.nan
        pf = pf.interpolate()

        # Vox
        vox = mmi - pf
        integral_vox = vox.sum()

        if oversmp and integral_vox < 500:
            continue

        # Barycenter
        index = vox.index
        unx_time = index.values.astype(np.int64) // 10 ** 9
        unx_sbc = (unx_time * vox).sum() / vox.sum()

        if unx_sbc < 0:  # TODO clean the situation
            print('Warning! unx_sbc < 0, in position:{0}'.format(coord))
            continue
        else:
            try:
                sbc = pd.to_datetime(unx_sbc * 10 ** 9)
            except (RuntimeError, Exception, ValueError):
                print('Warning! Datetime conversion went wrong, in position:{0}'.format(coord))
                continue

        # avoid unusual results
        if sbc.year not in range(pks.index[i].year - 1, pks.index[i + 1].year + 1):
            print('Warning! sbc not in a valid range, in position:{0}'.format(coord))
            continue
        # Season deviation standard
        sds = np.sqrt((np.square(unx_time) * vox).sum() / vox.sum() - np.square(unx_sbc))

        if np.isnan(sds):
            print('Warning! Season deviation standard is Nan, in position:{0}'.format(coord))
            continue

        sds_d = pd.Timedelta(sds, unit='s')

        # Update values
        row = pd.DataFrame([(pks.index[i], pks.index[i + 1], sbc, sds_d)],
                           columns=['sd', 'ed', 'sbc', 'ssd'],
                           index=[sbc])

        ts_table = ts_table.append(row)

    # Core
    try:
        for i in range(len(ts_table)):

            index = ts_table.iloc[i]['sbc']

            # specific mas
            mas = (ts_table.iloc[i]['ed'] - ts_table.iloc[i]['sd']) - (2 * prmts['mavmet'] * ts_table.iloc[i]['ssd'])

            if mas.days < 0:
                continue

            original = ps[ts_table.iloc[i]['sd']:ts_table.iloc[i]['ed']]

            timedelta = (ts_table.iloc[i]['ed'] - ts_table.iloc[i]['sd']) * 2/3

            buffered = ps[ts_table.iloc[i]['sd']-timedelta:ts_table.iloc[i]['ed']+timedelta]

            try:
                smth_crv = buffered.rolling(mas.days, win_type='boxcar', center=True).mean()
            except (RuntimeError, Exception, ValueError):
                print('Warning! Smoothed curv calculation went wrong, in position:{0}'.format(coord))
                err_table['sthcrv'] = 1
                return err_table

            smoothed = smth_crv[ts_table.iloc[i]['sd'] - timedelta:ts_table.iloc[i]['ed'] + timedelta]

            # baricenter for the smoothed one
            unx_time_smth = smoothed.index.values.astype(np.int64) // 10 ** 9
            unx_sbc_Y_smth = (unx_time_smth * smoothed).sum() / unx_time_smth.sum()

            # smoothed -= unx_sbc_Y_smth - unx_sbc_Y

            back = smoothed[:ts_table.iloc[i]['sbc']]\
                .shift(1, freq=pd.Timedelta(days=int(mas.days / 2)))[ts_table.iloc[i]['sd']:].dropna()

            forward = smoothed[ts_table.iloc[i]['sbc']:]\
                .shift(1, freq=pd.Timedelta(days=-int(mas.days / 2)))[:ts_table.iloc[i]['ed']].dropna()

            # original = original.reindex(pd.date_range(forward.index[0].date(), back.index[len(forward) - 1].date()))

            # region Plot
            # original.plot(style='r.')
            # back.plot(style='c:')
            # forward.plot(style='m:')
            # plt.show()
            # endregion

            sbd, sed, sbd_ts, sbd_ts = 4 * [None]

            # research the starting point of the season
            try:
                ge_sbd = original.ge(back)
                change_sbd = ge_sbd.rolling(window=2, min_periods=2).apply(lambda x: np.array_equal(x, [False, True]))
                sbd = change_sbd[change_sbd == 1][-1:]
                sbd_ts = pd.Series(original.loc[sbd.index], sbd.index)

            except (RuntimeError, Exception, ValueError):
                continue

            # research the end point of the season
            try:
                ge_sed = original.ge(forward)
                change_sed = ge_sed.rolling(window=2, min_periods=2).apply(lambda x: np.array_equal(x, [True, False]))
                sed = change_sed[change_sed == 1][-1:]
                sed_ts = pd.Series(original.loc[sed.index], sed.index)

            except (RuntimeError, Exception, ValueError):
                continue

            max_date = original.idxmax()

            if sed is None or sbd is None:
                continue
            else:
                # Season slope
                sslp = (sed_ts.values[0] - sbd_ts.values[0]) / (sed.index[0] - sbd.index[0]).days
                if not(abs(sslp) < 0.15):
                    continue

            try:
                # Start date of the season
                ts_table.set_value(index, 'sbd', sbd.index)

                # End date of the season
                ts_table.set_value(index, 'sed', sed.index)

                # Slope
                ts_table.set_value(index, 'sslp', sslp)

                # Season Lenght
                sl = (sed.index - sbd.index).to_pytimedelta()[0]
                ts_table.set_value(index, 'sl', sl)

                # Season permanet
                sp = sbd_ts.append(sed_ts).resample('D').asfreq().interpolate(method='linear')
                spi = sp.sum()
                ts_table.set_value(index, 'spi', spi)

                # Season Integral
                si = original.loc[sbd.index[0]:sed.index[0]].sum()
                ts_table.set_value(index, 'si', si)

                # Cyclic fraction
                cf = si - spi
                ts_table.set_value(index, 'cf', cf)

                # Active fraction
                af = original.loc[sbd.index[0]:max_date] - sp[:max_date]
                afi = af.sum()
                ts_table.set_value(index, 'afi', afi)

                # reference yr
                ref_yr = ts_table['sbd'].iloc[i]+((ts_table['sed'].iloc[i]-ts_table['sbd'].iloc[i])*2)/3
                ts_table.set_value(index, 'yr', ref_yr.year)

                # region Plotting
                # plt.ion()
                # original.plot(style='k:')
                # back.plot(style='c')
                # forward.plot(style='y')
                # sed_ts.plot(style='yo')
                # sbd_ts.plot(style='co')
                # plt.show()
                # endregion

            except (RuntimeError, Exception, ValueError):
                print('Error! populating ts_table went wrong, in year{1} @:{0}'.format(coord, index))
                continue

    except (RuntimeError, Exception, ValueError):
        print('Error! populating ts_table went wrong, in position:{0}'.format(coord))
        return

    # plt.ion()
    # # zoom the wresult window
    # mng = plt.get_current_fig_manager()
    # mng.window.state('zoomed')
    #
    # # plt.figure(1)
    # # ts_dek_unfix.plot(style='g--')
    # # ts_dek.plot(style='r--')
    #
    # # plt.figure(2)
    # # mng = plt.get_current_fig_manager()
    # # mng.window.state('zoomed')
    #
    # plt.subplot(511)
    # ts_dek_resc_unfix.plot(style='r--')
    # ts_clean.plot(style='r--')
    # pks.plot(style='r--')
    # pks.plot(style='ro')
    #
    # plt.subplot(512)
    # # ts_dek_resc_unfix.plot(style='k:')
    # ps.plot(style='g')
    # ps.loc[ts_table.loc[:, 'sbd'].dropna()].plot(style='co')
    # ps.loc[ts_table.loc[:, 'sed'].dropna()].plot(style='mo')
    #
    # # smoothed_crvs.plot(style='y.')
    # plt.subplot(513)
    # ts_table['sslp'].plot()
    #
    # plt.subplot(514)
    # ts_table['sl'].plot()
    #
    # plt.subplot(515)
    # ts_table['si'].plot()
    #
    # plt.show()
    # # plt.clf()

    return ts_table, err_table, season_lng
