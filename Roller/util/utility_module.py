import pandas as pd

def create_3D_linked_list(labels, numpy_array_3D, value_label):
    """returns a panel with interaction (x-axis) - value (y axis) - time (Z axis)"""
    windows_n = numpy_array_3D.shape[2]
    linked_list_3D ={}

    for i in xrange(windows_n):
        target_2D_array = numpy_array_3D[:,:,i]
        linked_list = create_linked_list(labels, target_2D_array, value_label)
        linked_list_3D[i] = linked_list
    return pd.Panel(linked_list_3D)

def create_linked_list(labels, numpy_array_2D, value_label):
    """labels and array should be in row-major order"""
    linked_list = pd.DataFrame({'regulator-target':labels, value_label:numpy_array_2D.flatten()})
    return linked_list
