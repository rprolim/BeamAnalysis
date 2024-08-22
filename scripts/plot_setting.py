def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    import matplotlib as mpl
    import numpy as np    
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def plot_customization(ax=None,params=None):
    import numpy as np    
    ax.set_xlim(params['xlims'][0],params['xlims'][1])
    ax.set_ylim(params['ylims'][0],params['ylims'][1])
    #ax.set_xscale(params['xscale']) <----it needs to include this one
    #ax.set_yscale(params['yscale']) <----it needs to include this one
    if not params['xlabel'][params['k']]==None: ax.set_xlabel(params['xlabel'][params['k']])
    if not params['ylabel'][params['k']]==None: ax.set_ylabel(params['ylabel'][params['k']])
    if params['disable_top_axis'   ]: ax.spines['top'   ].set_visible(False)    
    if params['disable_bottom_axis']: ax.spines['bottom'].set_visible(False)
    if params['disable_left_axis'  ]: ax.spines['left'  ].set_visible(False)
    if params['disable_right_axis' ]: ax.spines['right' ].set_visible(False)
    if params["disable_xtick"]: ax.set_xticks([])
    if params["disable_ytick"]: ax.set_yticks([])    
    
def plot_add_text(params=None, ax=None):
    import matplotlib.pyplot as plt
    font = {'family': params['plot_add_text_family'],
            'color':  params['plot_add_text_color'], 
            'weight': params['plot_add_text_weight'], 
            'size':   params['plot_add_text_size']}
    if ax!=None: ax.text(params['plot_add_text_x'], params['plot_add_text_y'], params['plot_add_text'], fontdict=font)
    else:       plt.text(params['plot_add_text_x'], params['plot_add_text_y'], params['plot_add_text'], fontdict=font)
    #horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
def colorbar_customization(params=None):
    import matplotlib.pyplot as plt    
    import numpy as np    
    params['cbar_ticks'] = np.linspace(params['fwhm'].min(),
                                       params['fwhm'].max(),
                                       params['cbar_ticks_length'])    
    if params['cbar_ticks_scale_invert']:params['cbar_ticks'] = np.flip(params['cbar_ticks'])
    params['cbar'] = plt.colorbar(params['cmap'],  
                                  location    = params['cbar_location'], 
                                  extend      = params['cbar_extend'], 
                                  pad         = params['cbar_pad'],
                                  shrink      = params['cbar_shrink'],
                                  aspect      = params['cbar_aspect'],
                                  format      = params['cbar_format'])
    params['cbar'].set_label(params['cbar_label'][params['k']],
                             rotation = params['cbar_label_rotation'], 
                             labelpad = params['cbar_label_labelpad'],
                             size     = params['cbar_label_size'],
                             weight   = params['cbar_label_weight'])
    params['cbar'].set_ticklabels(np.around(params['cbar_ticks'],decimals=params['cbar_ticks_decimals']))
    params['cbar'].ax.tick_params(labelsize = params['ticks_font_size'],
                                  rotation  = params['cbar_ticks_rotation'])
    
def colorbar_setting(params=None):
    import matplotlib as mpl
    import numpy as np        
    params['c']     = np.arange(1, params['n_lines'] + 1)
    params['ticks'] = np.linspace(params['c'].min(),
                                  params['c'].max(), 
                                  params['ticks_size'])    
    params['cbar_scale_norm'] = mpl.colors.Normalize(vmin=params['c'].min(), vmax=params['c'].max())
    params['cmap'] = mpl.cm.ScalarMappable(norm=params['cbar_scale_norm'], cmap=params['cmap'])
    params['cmap'].set_array([])    
    #return {"cmap":cmap, "c":c, "ticks":ticks}    

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    import matplotlib as mpl
    import numpy as np
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)
