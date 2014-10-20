#a thing that grabs different timepoints of data, can set window and step size.
import pandas as pd
class Roller:
    def __init__(self, file_path):
        """read the file, put it in a pandas dataframe? yeah that sounds fine"""
        self.raw_data = pd.read_csv(file_path,sep=" ")
        self.current_step = 0
        self.window_width = 3
        self.step_size = 1
        # get overall width of the time-course
        self.time_vec = self.raw_data['time'].unique()
        self.overall_width = len(self.time_vec)
        self.current_window = self.get_window()

    def get_window(self):
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

