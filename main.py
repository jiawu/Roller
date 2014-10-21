from Roller import Roller
from sklearn.preprocessing import Imputer

file_path = "compressed_katrina_data.txt"
gene_start_column = 5
roll_me = Roller(file_path, gene_start_column)

#get only TFs data, window size of 4
roll_me.set_window = 4

#get window of datapoints
#infer B coefficients

current_window = roll_me.get_window()

#impute missing values
imputer = Imputer(missing_values="NaN")
filled_matrix = imputer.fit_transform(current_window)
#each vector is a row
clf = linear_model.Lasso(alpha=0.2)

while (roll_me.next != "end"):



