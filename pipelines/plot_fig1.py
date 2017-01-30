import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#usage example: plot_fig1(tdr.window_list[0].explanatory_data)
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
  g_list = [0,1,2,3,4,5,6,7,8,9]

  fig, axes = plt.subplots(10,1, figsize=(30,20))
  # 10 plots

  for g in g_list:
    x_list = []
    y_list = []
    c_list = []
    n_list = [0,10,20,30,40]
    ni_list = [0,9,18,27,36]
    for n, ni in zip(n_list,ni_list):
      #for each gene (column, slice the rows)
      seg1 = data[n:n+5, g]
      x_list.append(seg1)
      y_list.append([x for x in range(ni,ni+5)])
      c_list.append(color_list[co])

      seg2 = data[n+4:n+9, g]
      x_list.append(seg2)
      y_list.append([x for x in range(ni+4,ni+9)])
      c_list.append(color_list[co+1])
      
      seg3 = data[n+8:n+10, g]
      x_list.append(seg3)
      y_list.append([x for x in range(ni+8,ni+10)])
      c_list.append(color_list[co+2])

      # repeat colors by adding 11
    for x,y,c in zip(x_list,y_list,c_list):
      axes[g].plot(y,x, color = c, lw=3)
    for k,v in axes[g].spines.items():
      v.set_visible(False)
      axes[g].set_xticks([])
      axes[g].set_yticks([])

  fig.savefig("fig1_test.pdf")


