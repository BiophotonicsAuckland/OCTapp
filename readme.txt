---------- Usage ----------
To start the program, simply double click OCTApp.py (if you have associated .py files with the python executable). Otherwise open a windows terminal in this folder (shift + right click -> open command window here) and type 'python OCTApp.py'.

Spectra files should be .txt files with the correct format.
Dispersion compensation, Linear and Non-linear filenames should be of the form '<filename>_st<num>_len<num>.txt'


---------- Troubleshooting ----------
If the application is doing something funky, first try closing the program, deleting 'settings.npy', and reopening the program.

If your computer stutters or locks up whilst doing 3D calculations, it could be one of two things:
a) Your computer has run out of memory (RAM) and starts to use your hard drive instead. This causes a lot of overhead and will drastically slow down calculation and responsiveness (but should still complete, given enough time).
You can verify that this is the case by opening Task Manager and looking for something like 'Available Physical Memory'.
	First, close any non-essential applications that might be using lots of memory. You can see what applications are using lots of memory by looking in Windows Task Manager
	Second, consider using the 'One En Face Image Only' setting instead. This uses substantially less memory, at the expense of having to recalculate if you want a different cross section.
	Thirdly, if possible upgrade your memory. 16GB Should be sufficient for the 3D calculations I was given to test the program with.

b) The program is usng too many threads for your CPU.
	By default the program uses n-1 threads, where n is the number of (logical) cores on the CPU. You can change this by opening the OCTApp.py file with a text editor and search for the line 'processCount = multiprocessing.cpu_count() - 1'.
	You can either hardcode the value by replacing the line with: 'processCount = 3' for three threads, or simply change the offset to '- 2', or '- 3', etc.

Sometimes running the program from within an IDE can make the 3D tab not work correctly (when calculation is started it stays stuck at 1%).
The problem has something to do with the type of python console the IDE uses. If possible, see if you can tell your IDE to use an external python console instead.
Otherwise this can be solved by running the program by simply double clicking the OCTApp.py file (if you have associated the .py file extension with python), or by opening a terminal in the folder and entering 'python OCTApp.py'
If all else fails, open the OCTApp.py file with a text editor, locate the line 'multithreadingEnabled = True' and change it to 'multithreadingEnabled = False'. This will be significantly slower (2-4x slower depending on number of CPU cores).



---------- Programming Info ----------
To build the standalone version, first (pip) install cx_freeze into the used python directory, if not already installed.
Once installed, open a terminal in this folder (shift + right click -> open command window here) and enter 'python setup.py build' into the windows terminal. 
This will make a folder called 'build'. Inside this folder will be an OCTApp.exe file. Everything in the build folder is required for the .exe to run.

Naming conventions:
	Functions and Variables:
		Functions and Variables are typically camelCase, but not all code was written by myself, so there may be sections that do not comply with this convention.
		Functions are typically prefixed with the name of the tab that they are used in (if used in multiple tabs, prefix is omitted).
		If a function is used in the default tab, the prefix is omitted.
		Standard Image Tab - <no prefix>
		3D Tab - threeD
		Image Processing Tab - processing
		
		Example: Image reset used in 3D tab:
		threeD_resetImage()
	
	UI Elements:
		UI elements are named as 'tab_name_type', where 'tab' is the name of the tab the element is in, and the type is what kind of UI element it is.
		The 'name' and 'type' components are all lowercase (not camelCase) to differentiate them from functions/variables.
		If an element is in the default tab ('standard image') the 'tab' part is omitted.
		
		A nonexhaustive list of UI element type names is given below:
		Button - button
		Checkbox - check
		Dropdown list - combo
		Slider - slider
		Progressbar - progressbar
		Radio checkbox - radio
		Spinbox (integer) - spin
		Spinbox (double-precision floating point) - dspin
		
		Example:
		A floating point spinbox in the 3D tab might be called:
		threeD_myspinbox_dspin