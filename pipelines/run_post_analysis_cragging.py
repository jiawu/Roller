from Swing_old.util.Analyzer import Analyzer
for i in range(1,6):
    path = '/projects/p20519/roller_output/optimizing_window_size/RandomForest/insilico_size10_' + str(i)

    #analyzer = Analyzer('/projects/p20519/roller_output/optimizing_window_size/RandomForest/ecoli')
    analyzer = Analyzer(path)
    outpath = '/projects/p20519/Swing/Swing/unittests/overall_df' + str(i) + '.csv'
    analyzer.overall_df.to_csv(outpath,index=False,sep=',')
