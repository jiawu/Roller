#a thing that grabs different timepoints of data, can set window and step size.
import pandas as pd
class Roller:
    def __init__(self, file_path,gene_start=None,gene_end=None):
        """read the file, put it in a pandas dataframe? yeah that sounds fine"""
        self.raw_data = pd.read_csv(file_path,sep=" ")
        self.current_step = 0
        self.window_width = 3
        self.step_size = 1
        # get overall width of the time-course
        self.time_vec = self.raw_data['time'].unique()
        self.overall_width = len(self.time_vec)
        if (gene_end != None):
            self.gene_end = gene_end
        else:
            self.gene_end = len(self.raw_data.columns)
        if (gene_start !=None):
            self.gene_start = gene_start
        else:
            self.gene_start = 0
        self.current_window = self.get_window()

    def get_window(self):
        raw_window = self.get_window_raw()
        only_genes=raw_window.iloc[:,self.gene_start:self.gene_end]
        return(only_genes)

    def get_window_raw(self):
        start_index = self.current_step
        end_index = start_index + self.window_width
        time_window = self.time_vec[start_index:end_index]
        data = self.raw_data[self.raw_data['time'].isin(time_window)]
        return data

    def next(self):
        end_index = self.current_step + self.window_width
        if (end_index <= self.overall_width) :
            self.current_step = self.current_step + self.step_size
            self.current_window = self.get_window()
            return(self.current_window)
        else:
            return("reached end")

    def set_window(self, width):
        self.window_width = width

    def set_step(self, step):
        self.step_size = step
    def reset(self):
        self.current_step = 0

