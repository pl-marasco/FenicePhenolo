import matplotlib.pyplot as plt


def plot(pxdrl):

    fig, axes = plt.subplots(3, 1, figsize=(20, 10))

    pxdrl.ts_raw.plot(ax=axes[0], style='k', title='Original Time Series for pixel in {0}'.format(pxdrl.position))

    pxdrl.ts_resc.plot(ax=axes[1], style='y', title='TS without masked values')

    # pxdrl.ts_cleaned.plot(ax=axes[1], style='y', title='TS interplated')

    gaps = pxdrl.ts_filtered[~(pxdrl.ts_filtered.shift(-1).notnull() & pxdrl.ts_filtered.shift(1).notnull())][1:-1]

    if len(gaps.values) > 0:
        gaps.plot(ax=axes[1], style='ro', title='TS without strong outlayers')

    pxdrl.ps.plot(ax=axes[2], style='g', title='TS smoothed with Savinsky Golet and braking points')

    if len(pxdrl.pks > 0):
        pxdrl.pks.plot(ax=axes[2], style='ro')

    plt.tight_layout()

    import math
    col = 3
    rows = math.ceil(len(pxdrl.phen) / col)

    fig, axes = plt.subplots(rows, col, figsize=(22, 10))
    plt.tight_layout()

    for i in range(0, len(pxdrl.phen)):
        phency = pxdrl.phen[i]
        plt.subplot(rows, col, i+1)

        if phency.buffered is not None:
            phency.mmc.plot(style='b', title='{}'.format(phency.ref_yr.values[0]))
        if phency.back is not None:
            phency.back.plot(style='-', color='gold')
        if phency.forward is not None:
            phency.forward.plot(style='-', color='olive')
        if phency.sbd is not None:
            phency.sbd.plot(style='r>',)
        if phency.sed is not None:
            phency.sed.plot(style='r<',)
        if phency.max is not None:
            phency.max.plot(style='rD')

    fig, axes = plt.subplots(5, 1, figsize=(22, 10))
    plt.tight_layout()

    plt.subplot(5, 1, 1)
    plt.tight_layout()
    if pxdrl.sl is not None:
        pxdrl.sl.plot(style='r', title='Season Lenght')
    plt.subplot(5, 1, 2)
    plt.tight_layout()
    if pxdrl.spi is not None:
        pxdrl.spi.plot(style='r', title='Season permanet')
    plt.subplot(5, 1, 3)
    plt.tight_layout()
    if pxdrl.si is not None:
        pxdrl.si.plot(style='r', title='Season Integral')
    plt.subplot(5, 1, 4)
    plt.tight_layout()
    if pxdrl.cf is not None:
        pxdrl.cf.plot(style='r', title='Cyclic fraction')
    plt.subplot(5, 1, 5)
    if pxdrl.afi is not None:
        pxdrl.afi.plot(style='r', title='Active fraction')

    plt.show()



