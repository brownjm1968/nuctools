import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as aa

__all__ = ['set_ticks_1axes','set_ticks_yt']

def set_ticks_1axes(axes_name):
	"""
	Used to set the ticks in the style Yaron Danon requires
	for his graduate students.
	
	Parameters
	----------
	axes_name : matplotlib axes object
	    This should be an axisartist axes object created like
	    in the example shown below
	
	Returns
	-------
	nothing
	
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> import mpl_toolkits.axisartist as aa
	>>> fig = figure(1)
	>>> ax1 = aa.Subplot(fig,111)
	>>> fig.add_subplot(ax1)
	>>> set_ticks_1axes(ax1)
	
	"""
	axes_name.axis["bottom","left"].major_ticks.set_tick_out(True)
	axes_name.axis["top","right"].major_ticks.set_tick_out(False)
	axes_name.axis["bottom","left"].minor_ticks.set_tick_out(True)
	axes_name.axis["top","right"].minor_ticks.set_tick_out(False)


def set_ticks_yt(axes_name1,axes_name2):
	"""
	Used to set the ticks in the style Yaron Danon requires
	for his graduate students. This in particular is designed
	for plotting yield and transmission on top of each other.
	
	Parameters
	----------
	axes_name1 : matplotlib axes object
	    This should be an axisartist axes object created like
	    in the example shown below
	axes_name2 : matplotlib axes object
	    This should be an axisartist axes object created like
	    in the example shown below
	
	axes_name : matplotlib axes object
	    This should be an axisartist axes object created like
	    in the example shown below
	
	Returns
	-------
	nothing
	
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> import mpl_toolkits.axisartist as aa
	>>> fig = figure(1)
	>>> ax1 = aa.Subplot(fig,211)
	>>> fig.add_subplot(ax1),fig.add_subplot(ax2)
	>>> set_ticks_yt(ax1,ax2)
	
	"""
	axes_name1.axis["bottom"].major_ticks.set_tick_out(False)
	axes_name1.axis["left"].major_ticks.set_tick_out(True)
	axes_name1.axis["bottom"].minor_ticks.set_tick_out(False)
	axes_name1.axis["left"].minor_ticks.set_tick_out(True)
	axes_name1.axis["top"].major_ticks.set_tick_out(False)
	axes_name1.axis["right"].major_ticks.set_tick_out(False)
	axes_name1.axis["top"].minor_ticks.set_tick_out(False)
	axes_name1.axis["bottom"].major_ticklabels.set_visible(False)
	#-----------------------------
	axes_name2.axis["bottom"].major_ticks.set_tick_out(True)
	axes_name2.axis["left"].major_ticks.set_tick_out(True)
	axes_name2.axis["bottom"].minor_ticks.set_tick_out(True)
	axes_name2.axis["left"].minor_ticks.set_tick_out(True)
	axes_name2.axis["top","right"].major_ticks.set_tick_out(False)
	axes_name2.axis["top","right"].minor_ticks.set_tick_out(False)
	axes_name2.axis["bottom"].major_ticklabels.set_visible(True)