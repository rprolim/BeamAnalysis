import numpy as np
def params_standard(
    beam_model = ['gaussian','cosine','jinc'],
    cmap_list = ['gist_heat','bone','Oranges'],
    xlim_min   = 0,
    xlim_max   = 4,
    ylim_min   = -50,
    ylim_max   = 1,
    ticks_size      = 9,
    ticks_font_size = 10,
    ticks_labels        = np.asarray("A B C D E F G H I J K L M N O P Q R S T U W V Y Z ".split(" ")),
    ticks_labels_true   = False,
    disable_xtick       = False,
    disable_ytick       = False,
    disable_top_axis    = True,
    disable_bottom_axis = False,
    disable_left_axis   = False,
    disable_right_axis  = True,
    cbar_ticks_decimals = 1,
    cbar_ticks_length   = 10,
    cbar_ticks_rotation = 45,
    cbar_orientation  = 'horizontal',
    cbar_location     = "top",
    cbar_extend       = 'both',
    cbar_pad          = 0.05,
    cbar_shrink       = 1.,
    cbar_aspect       = 70,
    cbar_format       = "%.1f",
    cbar_ticks_scale_invert=True,
    cbar_label_rotation = 0,
    cbar_label_labelpad = 20,
    cbar_label_size   = "small",
    cbar_label_weight = 'bold',
    cbar_label       = [r"$\theta_{\small \textrm{FWHM}}\ [\textrm{deg}]$",
                        r"$\theta_{\small \textrm{FWHM}}\ [\textrm{deg}]$",
                        r"$\theta_{\small \textrm{FWHM}}\ [\textrm{deg}]$"],
    xlabel = [r"$\theta \textrm{[deg]}$",
              r"$\theta \textrm{[deg]}$",
              r"$\theta \textrm{[deg]}$"],
    ylabel = [r'$\textrm{Beam}\ \textrm{pattern}\ \textrm{[dB]}$',
              None,
              None],
    plot_add_text_true   = False,
    plot_add_text        = r'$text$',
    plot_add_text_x      = 2,
    plot_add_text_y      = 0.65,
    plot_add_text_color  = "grey",
    plot_add_text_weight = 'normal',
    plot_add_text_family = 'serif',
    plot_add_text_size   = 16):  
    return {'beam_model':np.asarray(beam_model),
            "cmap_list":np.asarray(cmap_list),
            "ticks_size":ticks_size,  
            "ticks_labels":ticks_labels, 
            "ticks_labels_true":ticks_labels_true, 
            "ticks_font_size":ticks_font_size, 
            'cbar_ticks_length':cbar_ticks_length,
            'cbar_ticks_rotation':cbar_ticks_rotation,
            "cbar_orientation":cbar_orientation, 
            "cbar_location":cbar_location, 
            "cbar_extend":cbar_extend, 
            "cbar_pad":cbar_pad, 
            "cbar_shrink":cbar_shrink,
            "cbar_aspect":cbar_aspect, 
            "cbar_format":cbar_format,
            'cbar_ticks_decimals':cbar_ticks_decimals,
            "cbar_ticks_scale_invert":cbar_ticks_scale_invert,
            "cbar_label":np.asarray(cbar_label),
            'cbar_label_size'    :cbar_label_size,
            'cbar_label_weight'  :cbar_label_weight,
            'cbar_label_rotation':cbar_label_rotation,
            'cbar_label_labelpad':cbar_label_labelpad,
            "xlims": [xlim_min,xlim_max],
            "ylims": [ylim_min,ylim_max],
            "disable_xtick":disable_xtick,
            "disable_ytick":disable_ytick,
            "disable_top_axis":disable_top_axis,
            "disable_bottom_axis":disable_bottom_axis,
            "disable_left_axis":disable_left_axis,
            "disable_right_axis":disable_right_axis,
            'xlabel':np.asarray(xlabel),
            'ylabel':np.asarray(ylabel),
            "plot_add_text_true": plot_add_text_true,
            'plot_add_text':   plot_add_text,
            'plot_add_text_x': plot_add_text_x,
            'plot_add_text_y': plot_add_text_y,
            'plot_add_text_color': plot_add_text_color,
            'plot_add_text_weight': plot_add_text_weight,
            'plot_add_text_size': plot_add_text_size,
            'plot_add_text_family': plot_add_text_family            
         }    
