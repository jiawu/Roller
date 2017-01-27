import matplotlib.pyplot as plt

#data = tdr.window_list[0].explanatory_data
# 50 x 30

# 50 = 10* 5 perturbations
# 30 = 10 genes * 3 different lags (0, 1, 2)

#tdr.window_list[0].response_data
# 50 x 10

# 50 = 10 * 5 perturbations
# 10 genes at window 2

def plot_fig1(data):
  color_list = ['red', 'blue', 'green', 'gray', 'orange']
  # 5 colors

  # 4 timepoints a color. maximum is 10 timepoints.
  x_list = []
  y_list = []
  c_list = []
  co = 0

  n = 0

  #take each overall window, divide into segments based on overall_y
  #overall_y = [0,1,2,3,4,5,6,7,8,9]
  g = [0,1,2,3,4,5,6,7,8,9,10]

  fig, ax = plt.subplots(10,1)
  # 10 plots

  for g in g_list:
    n_list = [0,10,20,30,40]
    for n in n_list:
      #for each gene (column, slice the rows)
      seg1 = data[n:n+4, g]
      x_list.append(seg1)
      y_list.append([x for x in range(n,n+4)])
      c_list.append(color_list[co])

      seg2 = data[n+4:n+8, g]
      x_list.append(seg2)
      y_list.append([x for x in range(n+5,n+9)])
      c_list.append(color_list[co+1])

      seg3 = data[n+8:n+10, g]
      x_list.append(seg3)
      y_list.append([x for x in range(n+10,n+11)])
      c_list.append(color_list[co+2])

      # repeat colors by adding 11
    for x,y,c in zip(x_list,y_list,c_list):
      ax[g].plot(x,y, color = c)

  fig.savefig("fig1_test.png")


