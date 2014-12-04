import pandas as pd

import pdb
class Roller:
    """
    A thing that grabs different timepoints of data, can set window and step size.

    To do list:
        -i need to make two modules: a file processing module and a rolling module.
        -currently it only accepts tab separated files.
        -currently it looks for the "time" column. It looks for a hardcoded label"Time" which should be changed to a variable.
        -accept different kinds of files
    """
    def __init__(self, file_path, gene_start=None, gene_end=None, time_label="Time", separator = "\t"):
        """
        Initialize the roller object. Read the file and put it into a pandas dataframe
        :param file_path: file-like object or string
                        The file to read
        :param gene_start: int
        :param gene_end: int
        """
        # Read the raw data into a pandas dataframe object
        self.raw_data = pd.read_csv(file_path, sep=separator)
        self.raw_data = self.raw_data.dropna(axis=0, how='all')

        # Set roller defaults
        self.current_step = 0
        self.window_width = 3
        self.step_size = 1
        self.time_label = time_label
        #pdb.set_trace()
        # Get overall width of the time-course
        self.time_vec = self.raw_data[self.time_label].unique()
        self.overall_width = len(self.time_vec)
        if gene_end is not None:
            self.gene_end = gene_end
        else:
            self.gene_end = len(self.raw_data.columns)
        if gene_start is not None:
            self.gene_start = gene_start
        else:
            self.gene_start = 0
        self.current_window = self.get_window()

    def get_n_windows(self):
        total_windows = (self.overall_width - self.window_width)/(self.step_size)
        return total_windows

    def get_window(self):
        raw_window = self.get_window_raw()
        only_genes = raw_window.iloc[:, self.gene_start:self.gene_end]
        return only_genes

    def get_window_raw(self):
        start_index = self.current_step
        end_index = start_index + self.window_width
        time_window = self.time_vec[start_index:end_index]
        data = self.raw_data[self.raw_data[self.time_label].isin(time_window)]
        return data

    def next(self):
        end_index = self.current_step + self.window_width
        if end_index <= self.overall_width:
            self.current_step += self.step_size
            self.current_window = self.get_window()
            return self.current_window
        else:
            return "end"

    def set_window(self, width):
        self.window_width = width

    def set_step(self, step):
        self.step_size = step

    def reset(self):
        self.current_step = 0

    def fit(self):
        """
        Fit the rolling model

        :return: 3D matrix (n_windows, n_parents, n_children)
            Matrix of coefficients. Each slice is the beta coefficient matrix for a given window
        """
        file_path = "data/dream4/insilico_size10_1_timeseries.tsv"
        #file_path = "/Users/jjw036/Roller/goldbetter_model/goldbetter_data.txt"
        gene_start_column = 1
        time_label = "Time"
        separator = "\t"
        gene_end = None

        roll_me = Roller.Roller(file_path, gene_start_column, gene_end, time_label, separator)
        window_size = 5
        #get only TFs data, window size of 4
        roll_me.set_window(window_size)
        #impute missing values
        imputer = Imputer(missing_values="NaN")
        mpl.rcParams.update({'font.size':8})
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
        total_window_number = roll_me.get_n_windows()
        #todo: fix the below two lines, I don't need this to get the matrix size...
        current_window = roll_me.get_window()
        filled_matrix = imputer.fit_transform(current_window)
        n_genes = filled_matrix.shape[1]
        coeff_matrix_3d = np.empty((n_genes,n_genes,total_window_number))
        gene_names=list(current_window.columns.values)
        for nth_window in range(0, total_window_number):
        #loop gets the window, gets the coefficients for that window, then increments the window
        current_window = roll_me.get_window()

        #check if any values are completely blank
        sums = [current_window.iloc[:,x].sum() for x in range(0,10)]
        ind = np.where(np.isnan(sums))[0]
        current_window.iloc[:,ind]=0
        current_window = current_window *100
        filled_matrix = current_window.values
        #filled_matrix = imputer.fit_transform(current_window)
        current_lasso = LassoWrapper(filled_matrix)
        coeff_mat = current_lasso.get_coeffs(10)
        coeff_matrix_3d[:,:,nth_window] = coeff_mat
        #plot_figure(coeff_mat,nth_window,gene_names,gene_names,window_size)
        roll_me.next()