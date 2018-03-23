import numpy as np


def nodatafix(ts):
    # 252 cloud
    # 253 snow
    # 254 sea (not treated)
    # 251/255 invalids or data error

    # mean with nodata value included
    aggreg = ts.groupby([ts.index.month, ts.index.day]).median()

    # mean without nan included
    ntnlTs = ts.copy()
    ntnlTs[ntnlTs > 250] = np.nan
    ntnlAgg = ntnlTs.groupby([ts.index.month, ts.index.day]).median()
    ntnlAgg_min = ntnlTs.min()

    ntnlAggMth = None
    aggTriMth = None
    ntnlAggTriMth = None
    aggMth = None
    aggYR = None

    nanList = ts[ts > 250]

    for index, value in nanList.iteritems():
        aggVal = aggreg.loc[index.month, index.day]

        if value in [251, 252, 255]:
            if 250 < aggVal < 252.5 or aggVal >= 254.5:
                if aggMth is None:
                    aggMth = ts.groupby([ts.index.month]).median()
                if 250 < aggMth.loc[index.month] < 252.5 or aggMth.loc[index.month] >= 254.5:
                    if aggTriMth is None:
                        aggTriMth = ts.groupby([ts.index.quarter]).median()
                    if 250 < aggTriMth.loc[index.quarter] < 252.5:
                        if aggYR is None:
                            aggYR = ts.groupby([ts.index.year]).min()
                        # nanList.loc[index] = aggYR.loc[index.year]
                        nanList.loc[index] = aggYR.min()
                    else:
                        if ntnlAggTriMth is None:
                            ntnlAggTriMth = ntnlTs.groupby([ts.index.quarter]).median()
                        nanList.loc[index] = ntnlAggTriMth.loc[index.quarter]
                else:
                    if ntnlAggMth is None:
                        ntnlAggMth = ntnlTs.groupby([ts.index.month]).median()
                    nanList.loc[index] = ntnlAggMth.loc[index.month]
            elif 252.5 <= aggVal < 254.5:
                nanList.loc[index] = 0
            else:
                nanList.loc[index] = ntnlAgg.loc[index.month, index.day]

        elif value == 253:
            if aggVal > 250:
                nanList.loc[index] = ntnlAgg_min
            else:
                # nanList.loc[index] = ntnlAgg.loc[index.month, index.day]
                nanList.loc[index] = np.NaN

    ts.update(nanList)

    ts[ts == 253] = np.NaN
    ts = ts.interpolate()

    return ts
