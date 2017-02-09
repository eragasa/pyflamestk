import pyflamestk.pareto as pareto
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

class ParetoPostProcessor:
    def __init__(self,sim_results):
        assert type(sim_results),pareto.SimulationResults
        self._sim_results = sim_results
    
    @property
    def sim_results(self): return self._sim_results

    @sim_results.setter
    def sim_results(self,sim_results):
        assert type(sim_results),pareto.SimulationResults
        self._sim_results = sim_results
        
    @property
    def free_parameters(self): return self._free_params

    @free_parameters.setter
    def free_parameters(self, param_dict):
        self._free_params = param_dict

def simresults_create_2d_pareto_plot(simresults,
                          qoi_name_1,
                          qoi_name_2,
                          error_type = "",
                          show_dominated = True,
                          show_pareto = True, 
                          show_culled = True,
                          show_2d_pareto_curve = True,
                          show_2d_culled_curve = True):   
        """
        Arguments:
        
        qoi_name_1 (string) - the name of the first quantity of interest to be
            plotted on the x-axis
        qoi_name_2 (string) - the name of the second quantity of interest to be
            plotted on the y-axis
        error_type (string: "" ) - not implemented yet
        show_dominated (bool, True) - shows dominated points when set to true
        show_pareto (bool,True) - shows Pareto points when set to true.
        show_culled (bool,True) - shows the remaining Pareto points when the 
            Pareto set is culled by by performance requirements.
        show_2d_pareto_curve (bool,True) - shows a 2d slice of the Pareto curve
            when set to True.
        show_2d_pareto_cuve (bool,True) - shows a 2d slice of the culled Pareto
            curve when set to false.
            
        Returns:
        
        Nothing
            
        Notes:
        
        When the dominated points are plotted, it actually plots all the points
        in the dataset.  However, these points are then plotted over by the
        Pareto points.
        """
        
        
        x_label = qoi_name_1
        y_label = qoi_name_2
                                   
        x_data_all = simresults.get_data(x_label, 'all')
        y_data_all = simresults.get_data(y_label, 'all')
        
        if show_pareto == True:
            x_data_pareto = simresults.get_data(x_label, 'pareto')
            y_data_pareto = simresults.get_data(y_label, 'pareto')

            if show_2d_pareto_curve == True:
                pareto_2d = pareto.pareto_frontier_2d(x_data_all,
                                                      y_data_all,
                                                      maxX = False, 
                                                      maxY= False) 
            
        if show_culled == True:
            x_data_cull = simresults.get_data_by_name(x_label, 'pareto_cull')
            y_data_cull = simresults.get_data_by_name(y_label, 'pareto_cull')

            if show_2d_pareto_curve == True:
                cull_2d = pareto.pareto_frontier_2d(x_data_cull,
                                                    y_data_cull,
                                                    maxX = False, 
                                                    maxY= False)         
        # plot results
        plt.figure()
    
        if show_dominated == True:
            c = 'b'
            plt.scatter(x_data_all, y_data_all, color = c)
            
        if show_pareto == True:
            c = 'y'
            plt.scatter(x_data_pareto,y_data_pareto, color=c)
            if show_2d_pareto_curve == True:
                plt.plot(pareto_2d[0],pareto_2d[1],color=c)

        if show_culled == True:
            c = 'g'
            plt.scatter(x_data_cull,  y_data_cull, color = c)
            if show_2d_pareto_curve == True:
                plt.plot(cull_2d[0],cull_2d[1],color=c)

        # determine axis
        xmin = min(pareto_2d[0])
        xmax = max(pareto_2d[0])
        ymin = min(pareto_2d[1])
        ymax = max(pareto_2d[1])
        plt.axis([xmin,xmax,ymin,ymax])
        
        # add labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        
        # display graph
        plt.show()


def simresults_create_histogram(simresults,
                                name, 
                                ds_type,
                                err_type = 'abs_err'):
    # check type
    assert type(simresults), pareto.SimulationResults
    assert type(name), str
    assert type(ds_type), str
    assert type(err_type), str

    # get data
    data = simresults.get_data(name,ds_type,err_type)
    n_data = data.shape[0]

    # remove outliers
    outliers = is_outlier(data)
    i_outliers = [i for i in range(n_data) if outliers[i] == False]
    data = data[i_outliers]

    m_data = data.shape[0]
    # plot figure
    print('\t mdata/ndata={}/{}'.format(m_data,n_data))
    plt.figure()
    plt.hist(data)
    plt.show()

def simresults_create_param_histograms(simresults,ds_type):
    # check type
    assert type(simresults), pareto.SimulationResults
    assert type(ds_type), str

    # get parameter names
    p_names = simresults.parameter_names

    for p in p_names:
        print('param_histogram:',p,'.',ds_type)
        simresults_create_histogram(simresults,p,ds_type)

def simresults_create_qoi_histograms(simresults,ds_type):
    # check type
    assert type(simresults), pareto.SimulationResults
    assert type(ds_type), str

    # get parameter names
    q_names = simresults.qoi_names

    for q in q_names:
        print('qoi_histogram',q,'.',ds_type)
        simresults_create_histogram(simresults,q,ds_type)

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 

        original code: http://stackoverflow.com/questions/11882393/matplotlib-disregard-outliers-when-plotting
        date: 10/19/2016
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh