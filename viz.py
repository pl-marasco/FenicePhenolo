import matplotlib.pyplot as plt


def plot(pxdrl):

    fig, axes = plt.subplots(3, 1)

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

    plt.show()
