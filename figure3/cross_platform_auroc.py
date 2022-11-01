import argparse
import os

from collections import Counter

from matplotlib import pyplot as plt
import seaborn as sns

col_types = [
    # hard-coded datatypes of the columns in the input file
    # sorry bout that :(
    int, str, str, int, str, str,
    int, float, str, str, str, str,
    list, float, float,
    str, str, str, str, str,
    int, # last column is DP_ALT and will be calculated on the fly
    str # the original groups string
]

group_defaults = {
    'ILLUMINA': [
        'BEI',
        'Galaxy',
        'Gilead',
        'LANL',
        'Lilly',
        'NCBI',
        'VIR',
    ],
    'OXFORD_NANOPORE': [
        'BEI',
        'Galaxy',
        'Gilead',
        'LANL',
        'NCBI',
    ]
}

tp_criteria = {
    'ILLUMINA': 4,
    'OXFORD_NANOPORE': 3
}

param_ranges = {
    # CAUTION: currently, all ranges need to start at 0 or false-positive
    # rates will be normalized incorrectly when -n/--normalize-fp
    # is in use
    'DP_ALT': [0, 10, 20, 30, 50, 100, 200, 400, 1000, 2000, 4000],
    'DP': [0, 10, 20, 30, 50, 100, 200, 400, 1000, 2000, 4000],
    'AF': [n/100 for n in range(0, 101, 10)]
}

# Galaxy has only analyzed a subset of the ActivTrace data.
# We read the analyzed accession numbers from an accompanying file.
with open(
    os.path.join(os.path.dirname(__file__), 'Galaxy_analyzed.txt')
) as gx_accs:
    gx_analyzed = {acc.strip() for acc in gx_accs}


def rocplot(rates, ofile, kind, title):
    # sns.set_theme(rc={"figure.figsize":(16, 8)}) # set figure width and height
    sns.set_context('paper')
    sns.set_style('whitegrid')
    x = []
    y = []
    labels = []
    hues = []
    for group, data in rates.items():
        rates = sorted(data.items(), reverse=True)
        # use param_values as labels
        labels.extend([r[0] for r in rates])
        # use true-positive rates as y
        y.extend([r[1][0] for r in rates])
        # false-positive rates as x
        x.extend([r[1][1] for r in rates])
        hues.extend([group]*len(rates))
    if kind == 'line':
        kwargs = {'estimator': None}
    else:
        kwargs = {}
    g = sns.relplot(
        x=x, y=y,
        kind=kind, hue=hues, marker='o', palette='colorblind',
        height=8, aspect=1.5, **kwargs
    )
    plt.xlabel('Discordant calls rate')
    plt.ylabel('Concordant calls rate')
    g._legend = None
    g.axes[0, 0].set_xlim(0, 1.05)
    g.axes[0, 0].set_ylim(0, 1.02)
    g.axes[0, 0].set_xticks(ticks=[0.2, 0.4, 0.6, 0.8, 1.0])
    g.axes[0, 0].set_yticks(ticks=[0.2, 0.4, 0.6, 0.8, 1.0])
    g.axes[0, 0].legend(loc='lower right', fontsize='large', title=title)
    old_coords = (0, 0)
    # put labels on some of our data points
    # for DP and DP_ALT label multiples of 50
    # for AF label mutiples of 0.1
    for param_value, xpos, ypos in zip(labels, x, y):
        if (xpos > 0.99 and ypos > 0.99) or (xpos < 0.01 and ypos < 0.01):
            continue
        if param_value and (
            not param_value % 50
            or (not (param_value * 100) % 10)
        ):
            if (xpos, ypos) != old_coords:
                plt.annotate(xy=(xpos + 0.002, ypos - 0.015), text=str(param_value), fontsize=10)
                old_coords = (xpos, ypos)
    plt.savefig(ofile)


parser = argparse.ArgumentParser()
parser.add_argument(
    'ifile',
    help='Combined results csv file to use as input'
)
parser.add_argument(
    'ofile',
    help='Save the resulting ROC plot here'
)
parser.add_argument(
    '--platform', '--pl', choices=['ILLUMINA', 'OXFORD_NANOPORE'],
    required=True,
    help='Plot data for this platform'
)
parser.add_argument(
    '-p', '--parameter', choices=['AF', 'DP', 'DP_ALT'], default='DP_ALT',
    help='Generate ROC curve for this parameter (default: DP_ALT)'
)
parser.add_argument(
    '-t', '--tp-criterion', type=int, default=0,
    help='In -c groups mode, consider a variant found by this many groups a '
         'true-positive call (default: 4 for ILLUMINA, 3 for OXOFORD_NANOPORE).'
)
parser.add_argument(
    '-f', '--fp-criterion', type=int, default=0,
    help='In -c groups mode, apply this number as a separate threshold for '
         'false-positive calls (default: consider any non-true-positive a '
         'false-positive call).'
)
parser.add_argument(
    '-n', '--normalize-fp', action='store_true',
    help='Normalize false-positive rates (make them reach 1 for parameter '
         'value 0). The non-normalized default is the false-positive rate '
         'with regard to the total pool of likely false-positives across '
         'all groups.'
)
parser.add_argument(
    '--normalize-tp', action='store_true',
    help='Normalize true-positive rates (make them reach 1 for parameter '
         'value 0). The non-normalized default is the true-positive rate '
         'with regard to the total pool of likely true-positives across '
         'all groups.'
)
parser.add_argument(
    '-v', '--vartype', choices=['SNP', 'InDel'],
    help='Generate ROC curve for this type of variants only'
)
parser.add_argument(
    '-k', '--plot-type', choices=['scatter', 'line'], default='scatter',
    help='The kind of plot (scatter or lines plot) to generate'
)
parser.add_argument(
    '-c', '--classifier', choices=['groups', 'platforms'], default='groups',
    help='How to classify calls as true- or false-positive. If platforms is '
         'specified, a true-positive call is one that is made on both '
         'platforms, otherwise true-positives are defined as calls made by '
         'at least --tp-criterion groups.'
)

args = parser.parse_args()
tp_criterion = args.tp_criterion or tp_criteria[args.platform]
if args.classifier == 'platforms':
    groups = group_defaults['ILLUMINA']
else:
    groups = group_defaults[args.platform]

with open(args.ifile) as i:
    hdr_cols = i.readline().strip().split(',')
    data = []

    for line in i:
        fields = line.strip().split(',')
        if fields[hdr_cols.index('Group')] not in groups:
            # This record describes a call by a group we are not interested in
            continue
        # add DP_ALT on the fly
        fields.append(
            round(
                int(fields[hdr_cols.index('DP')])
                * float(fields[hdr_cols.index('AF')])
            )
        )
        # split Groups column and keep only groups of interest
        # store original Groups string as the last field
        fields.append(fields[hdr_cols.index('Groups')])
        fields[hdr_cols.index('Groups')] = [
            g for g in fields[hdr_cols.index('Groups')].split('_')
            if g in groups
        ]
        hdr_cols.append('DP_ALT')
        assert len(fields) == len(col_types), '{}: {}'.format(len(fields), col_types)
        data.append([t(f) for t, f in zip(col_types, fields)])
param_col = hdr_cols.index(args.parameter)
data.sort(key=lambda x: x[param_col])

if args.classifier == 'groups':
    tps = set()
    tps_gx = set()
    fps = set()
    for record in data:
        if record[hdr_cols.index('platform')] != args.platform:
            continue
        if args.vartype and record[hdr_cols.index('type')] != args.vartype:
            continue
        if args.fp_criterion and len(
            record[hdr_cols.index('Groups')]
        ) <= args.fp_criterion:
            fps.add(record[hdr_cols.index('var')])
            continue
        if args.platform == 'OXFORD_NANOPORE' or record[hdr_cols.index('Acc')] in gx_analyzed:
            if len(record[hdr_cols.index('Groups')]) >= tp_criterion:
                tps.add(record[hdr_cols.index('var')])
                tps_gx.add(record[hdr_cols.index('var')])
            elif not args.fp_criterion:
                fps.add(record[hdr_cols.index('var')])
        else:
            # This sample has not been analyzed by Galaxy and Lilly.
            # With only 5 groups instead of 7, we relax the tp_criterion by one.
            if len(record[hdr_cols.index('Groups')]) >= tp_criterion - 1:
                tps.add(record[hdr_cols.index('var')])
            elif not args.fp_criterion:
                fps.add(record[hdr_cols.index('var')])
    fp_all = fp_all_gx = len(fps)
    tp_all = len(tps)
    tp_all_gx = len(tps_gx)

data_by_platform_group_sample = {
    'ILLUMINA': {group: {} for group in groups},
    'OXFORD_NANOPORE': {group: {} for group in groups}
}
platform_col = hdr_cols.index('platform')
group_col = hdr_cols.index('Group')
sample_col = hdr_cols.index('biosample')
var_col = hdr_cols.index('var')
key_cols = sorted([platform_col, group_col, sample_col], reverse=True)
for record in data:
    sample_id = record[sample_col]
    var = (record[3], record[4], record[5])
    if sample_id in data_by_platform_group_sample[record[platform_col]][
        record[group_col]
    ]:
        data_by_platform_group_sample[record[platform_col]][record[group_col]][
            sample_id
        ][var] = record
    else:
        data_by_platform_group_sample[record[platform_col]][record[group_col]][
            sample_id
        ] = {var: record}
    # no need to retain record columns that have become keys in our dict
    for col in key_cols:
        try:
            del record[col]
        except TypeError:
            print(record, col, record[col])
            raise
# recalculate the record param_col now that we've shrunk records
for col in key_cols:
    del hdr_cols[col]
param_col = hdr_cols.index(args.parameter)
rates = {}

if args.classifier == 'platforms':
    samples = {
        'ILLUMINA': {
            sample
            for group_data in
            data_by_platform_group_sample['ILLUMINA'].values()
            for sample in group_data
        },
        'OXFORD_NANOPORE': {
            sample
            for group_data in
            data_by_platform_group_sample['OXFORD_NANOPORE'].values()
            for sample in group_data
        }
    }

    sample_vars = {
        'ILLUMINA': {
            (sample, var)
            for group_data in
            data_by_platform_group_sample['ILLUMINA'].values()
            for sample, sample_data in group_data.items()
            for var, record in sample_data.items()
            if (
                not args.vartype
            ) or (
                record[hdr_cols.index('type')] == args.vartype
            )
        },
        'OXFORD_NANOPORE': {
            (sample, var)
            for group_data in
            data_by_platform_group_sample['OXFORD_NANOPORE'].values()
            for sample, sample_data in group_data.items()
            for var, record in sample_data.items()
            if (
                not args.vartype
            ) or (
                record[hdr_cols.index('type')] == args.vartype
            )
        }
    }

    acc_col = hdr_cols.index('Acc')
    sample_vars_gx = {
        'ILLUMINA': {
            (sample, var)
            for group_data in
            data_by_platform_group_sample['ILLUMINA'].values()
            for sample, sample_data in group_data.items()
            for var, record in sample_data.items()
            if record[acc_col] in gx_analyzed and (
                not args.vartype
                or record[hdr_cols.index('type')] == args.vartype
            )
        },
        'OXFORD_NANOPORE': {
            (sample, var)
            for group_data in
            data_by_platform_group_sample['OXFORD_NANOPORE'].values()
            for sample, sample_data in group_data.items()
            for var, record in sample_data.items()
            if record[acc_col] in gx_analyzed and (
                not args.vartype
                or record[hdr_cols.index('type')] == args.vartype
            )
        }
    }

    if args.platform == 'ILLUMINA':
        other_platform = 'OXFORD_NANOPORE'
    else:
        other_platform = 'ILLUMINA'

    tps = set()
    tps_gx = set()
    fps = set()
    fps_gx = set()
    
    for group, group_data in data_by_platform_group_sample[
        args.platform
    ].items():
        for sample, sample_data in group_data.items():
            if sample not in samples[other_platform]:
                # this sample has not been analyzed on both platforms
                continue
            for var, var_data in sample_data.items():
                var_type = var_data[hdr_cols.index('type')]
                if args.vartype and var_type != args.vartype:
                    continue
                if (sample, var) in sample_vars[other_platform]:
                    # this variant has been called by someone also on the
                    # other platform => lets assume it's real
                    tps.add((sample, var))
                    if (
                        (sample, var) in sample_vars_gx[args.platform]
                    ) and (
                        (sample, var) in sample_vars_gx[other_platform]
                    ):
                        # Galaxy has also analyzed the sample
                        tps_gx.add((sample, var))
                else:
                    # nobody has made that call on the other platform.
                    if group == 'Galaxy' or group == 'Lilly':
                        if (sample, var) not in sample_vars_gx[other_platform]:
                            # but this is a call that Galaxy made without
                            # having seen the data from the other platform.
                            # Unless the same call appears for another group
                            # do not consider it a false-positive (could be
                            # a Galaxy-specific but cross-platform call).
                            continue
                    if group == 'NCBI' and var_type == 'InDel':
                        # NCBI doesn't call indels from ONT data so
                        # unless the same call appears for another group
                        # do not consider it a false-positive (could be
                        # a NCBI-specific but cross-platform call).
                        assert args.platform == 'ILLUMINA'
                        continue
                    if group not in group_defaults[other_platform]:
                        # This group hasn't called any variants from ONT data.
                        # Unless the same call appears for another group
                        # do not consider it a false-positive (could be
                        # a group-specific but cross-platform call).
                        continue
                    if (sample, var) in sample_vars_gx[args.platform]:
                        # This is a false-postive call on the platform,
                        # *and* Galaxy had a chance to call it cause it
                        # did analyze the sample causing it.
                        fps_gx.add((sample, var))
                    fps.add((sample, var))
    tp_all = len(tps)
    tp_all_gx = len(tps_gx)
    fp_all = len(fps)
    fp_all_gx = len(fps_gx)
    print('Counts:', tp_all, fp_all)


param_col = hdr_cols.index(args.parameter)
for platform, platform_data in data_by_platform_group_sample.items():
    if platform != args.platform:
        continue
    for group, group_data in platform_data.items():
        param_min = 10**99
        #    if group == 'Gilead':
        #        continue
        counts = {}
        for threshold in param_ranges[args.parameter]:
            true_counts = 0
            false_counts = 0
            for sample, sample_data in group_data.items():
                for var, record in sample_data.items():
                    if record[param_col] >= threshold and (
                        not args.vartype
                        or record[hdr_cols.index('type')] == args.vartype
                    ):
                        if record[param_col] < param_min:
                            param_min = record[param_col]
                        if args.classifier == 'groups':
                            if group == 'Galaxy' or group == 'Lilly':
                                if record[hdr_cols.index('var')] in tps_gx:
                                    true_counts += 1
                            else:
                                if record[hdr_cols.index('var')] in tps:
                                    true_counts += 1
                            if record[hdr_cols.index('var')] in fps:
                                false_counts += 1
                        else:
                            if group == 'Galaxy' or group == 'Lilly':
                                if (sample, var) in fps_gx:
                                    false_counts += 1
                                if (sample, var) in tps_gx:
                                    true_counts += 1
                            else:
                                if (sample, var) in fps:
                                    false_counts += 1
                                if (sample, var) in tps:
                                    true_counts += 1
            counts[threshold] = [
                true_counts,
                false_counts,
            ]
        if group == 'Galaxy' or group == 'Lilly':
            tp_total = tp_all_gx
            fp_total = fp_all_gx
        else:
            tp_total = tp_all
            fp_total = fp_all
        if args.normalize_fp:
            fp_div = counts[next(iter(param_ranges[args.parameter]))][1]
            if fp_div == 0:
                print('No false-positives for the data from', group)
                # exclude it from plotting
                continue
        else:
            fp_div = fp_total
        if args.normalize_tp:
            tp_div = counts[next(iter(param_ranges[args.parameter]))][0]
            if tp_div == 0:
                print('No true-positives for the data from', group)
                # exclude it from plotting
                continue
        else:
            tp_div = tp_total
        rates[group] = {
            t: [c[0]/tp_div, c[1]/fp_div] for t, c in counts.items()
        }
        print(group)
        print(tp_total, fp_total, param_min)

rocplot(
    rates,
    args.ofile,
    args.plot_type,
    'Effect of {0} on {1}'.format(args.parameter, args.platform)
)
