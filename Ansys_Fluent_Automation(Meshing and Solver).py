# Imported Packages
import ansys.fluent.core as pyfluent
import math
import time

################################################################################################################################################################################
# Input Section
################################################################################################################################################################################
#######################################################
# Inputs for the working paths of the meshing/Setup/Data files to contain in the loop.
#######################################################
# for this section Make sure the slashes are / not \
Parent_Folder = "C:/File_Path_Here"
# Note, all sub folders and cad files contained in this parrent folder must follow the following convention for automation:
# a<Value For a>b<Value For b>l<Value for l>IDAM<Value for IDAM>IDAI<Value for IDAI>AMNUM<Value for AMNUM>FC<value of the combustor fillet radius (top edge)>EL<Value of the exaust length in mm>
# An example of this is: "a40b25l100IDAM6IDAI18AMNUM8FC10EL30"
List_File_Names_To_Run = ["a20b13l100IDAM4IDAI18AMNUM8FC10EL30"] # In formate ["FileName1", "FileName2", etc]

#######################################################
# Meshing Setup Parameters
#######################################################
Meshing_Run = False
Meshing_Show_GUI = True

# Table of standard mesh sizes that will be used when generating the meshes
Table_Of_standard_Mesh_Sizes_Tube_Radius = [1,5,10,15,20,30,40,50,60,70,80]
Default_Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides = [0.05,0.15,0.35,0.45,0.55,0.70,0.85,0.95,1,1.05,1.1]
Meshing_Element_Size_Increase_Ammount = 0.025
Meshing_Element_Size_Decrease_Ammount = 0.025
# The above parameter is the increase in the element size that will take place in the while loop untill the mesh is acceptable Helps reduce element count in full model

# Boundry layer setup values
Boundry_Layer_Settings_First_Height = 0.2 # Given in mm
Boundry_Layer_Settings_Number_Of_layers = 5 # Kept consitant in making table
Boundry_Layer_Settings_Rate = 1.2 # kept consistant in making table
# Note the above variables can be adjusted with further reasurch to make better meshes for the geometries

# Periodic Boundry condition faces mesh size refinement
Meshing_Periodic_Boundry_Element_Size = 0.4 # given in mm
Meshing_Combustor_Outlet_Face_Element_Size = 0.4 # given in mm

#######################################################
# Solver Setup Parameters
#######################################################
Solver_Run = True
Solver_Run_Combustion = False
Solver_GPU_Enhanced = True
# For modeling combustion in 2023 R2 Ansys cannot compute combustion on a gpu.
# Therefore, Solver_GPU_Enhancd will auto disable if combustion is set to True.
Solver_Show_GUI = True
Solver_Session_Coupled_Solver_Time_Scale_Factor_Global = 1 # Default is 1 but small numbers may lead to better convergence.
Solver_Session_Coupled_Solver_Time_Scale_Factor_Energy = 0.1 # Default is 1
Solver_Session_Coupled_Solver_Time_Scale_Factor_Mean_Mixture_Fraction = 0.1 # Default is 1
Solver_Max_Iteration_Count = 3000

#######################################################
# Solver pathlines report options
#######################################################
Solver_Pathlines_Faces_Output = ["air_inlet_face","ammonia_inlet_faces"] # The faces that pathlines will start from
# Note the above names have to be exact. Use the following list as default: ["air_inlet_face","ammonia_inlet_faces"]
Solver_Pathlines_Steps = 2000 # steps for all path lines note default of 1000

#######################################################
# Inputs of fuel properties (Ammonia (gas) for this project) (At 1 Atm at 300K)
#######################################################
Fuel_Total_Flow_Rate = 0.000672043010752688 # Given in Kg per s
Fuel_Denisty = 0.6898  # Given in Kg per meter cubed
Fuel_Dynamic_Viscosity = 0.00001 # Given in Pascal Secounds
Fuel_Inlet_Tempurature = 300 # Given in Kelvin

#######################################################
# Inputs of air properties (Properties at 1Atm at 300K)
#######################################################
Air_Total_Flow_Rate = 0.00338825408963781 # Given in Kg per s
Air_Density = 1.177 # Given in Kg per meter cubed
Air_Dynamic_Viscosity = 0.00001845 # Given in Pascal Secounds
Air_Inlet_Tempurature = 300 # Given in Kelvin

#######################################################
# Combustion Setup Parameters
#######################################################
Combustion_Operating_Pressure = 101325
# Note the folder described bellow must be in the Parent folder described for the script to work
Combustion_Kinetics_Thermodynamics_Tranport_Folder_Name = "Folder_Name_Here"
# The file names of the chemical property files
Combustion_File_Name_Kinetics = "File_Name_Here.inp"
Combustion_File_Name_Thermodynamics = "File_Name_Here.dat"
Combustion_File_Name_Transport = "File_Name_Here.dat"

# The selected combustion species
# Common valid combustion species are:
# Methane - "ch4"
# Hydrogen - "h2"
# Ammonia - "nh3"
# As long as the species fuel is in the Chemkin files above that name will work
# The above names are just the ones of interest
Combustion_Fuel_Species_Name = "nh3"
# Note that this script assumed one fuel species burning in are where only o2 and n2 are considered as part of the mixture
Combustion_Equivalance_Ratio = 1.2                                  

# NOx parameters
Combustion_NOx_Enable = False
Combustion_NOx_O_Rad = 1 # Model Choice for O radicals (1 oxygen atom): 0 = equilibrium, 1 = partial-equilibrium, 2 = instantaneous
Combustion_NOx_OH_Rad = 0 # Model Choice for OH radicals (1 oxygen and hydrogen atom): 0 = equilibrium, 1 = partial-equilibrium, 2 = instantaneous

# Species Report settings
# displays the species mass weighted average at the outlet of the combustor
# This also enables a convergence criteria for all the species.
Combustion_Report_Refinement_Value = 0.001
Combustion_Report_Fuel = True # Inlcudes fuel species specied in the variable Combustion_Fuel_Species_Name
Combustion_Report_NOx = True # Includes no for hydrogen and methane. Includes no, no2 and nh2 for ammonia.
Combustion_Report_Carbons = False # Includes co2, co and c
# To add other species to report add in the list bellow. In the formate "species1","species2",etc
Combustion_Report_Species = [] # Valid species inputs are bellow this line:
# Note lower case letters
    # For Methane and Hydrogen Kinetics and Thermodynamics Data:
    #h2      h       o       o2      oh      h2o     ho2    
    #c       ch      ch2     ch3     ch4     co      co2     
    #ch2o    c2h2    c2h4    c2h6
    #nh3     no      hcn     n2

    # For Ammonia Kinetics and Thermodynamics Data:
    #no     nh3
    #h2     o2       h       o       oh      ho2     h2o    h2o2
    #nh2    nh       n
    #h2no   hnOh     hno     no2     hono    n2o
    #n2h4   n2h3     n2h2    h2nn 
    #n2

# Tempurature contour plots and convergence criteria settings
Combustion_Contour_Plot_Tempurature = True # Plots the 2d tempurature contour of the combustor
Combustion_Convergence_Tempurature = 0.001 # Convergence criteria value for tempurature at the outlet

# Contour plots settings
# Species specific contour settings
Combustion_Contour_Plot_Species = True # Enables plotting of 2d contor of species
Combustion_Contour_Plot_Fuel = True # Inlcudes fuel species specied in the variable Combustion_Fuel_Species_Name
Combustion_Contour_Plot_NOx = True # Includes no for hydrogen and methane. Includes no, no2 and nh2 for ammonia.
Combustion_Contour_Plot_Carbons = False # Includes co2, co and c
# To add other species to report add in the list bellow. In the formate "species1","species2",etc
Combustion_Contour_Species_To_Plot = ["n2","o2"] # Valid species inputs are bellow this line:
# Note lower case letters
    # For Methane and Hydrogen Kinetics and Thermodynamics Data:
    #h2      h       o       o2      oh      h2o     ho2    
    #c       ch      ch2     ch3     ch4     co      co2     
    #ch2o    c2h2    c2h4    c2h6
    #nh3     no      hcn     n2

    # For Ammonia Kinetics and Thermodynamics Data:
    #no     nh3
    #h2     o2       h       o       oh      ho2     h2o    h2o2
    #nh2    nh       n
    #h2no   hnOh     hno     no2     hono    n2o
    #n2h4   n2h3     n2h2    h2nn 
    #n2

# The following line enables monitoring points in line with the outlet, inlet and half the inlet radius
Combustion_Monitor_Points_Enable = True # enable or not
Combustion_Monitor_Points_Convergence_Number = 0.001 # Convergence criteria number
Combustion_Monitor_Points_Ammount_In_Line = 2 # Ammount of points to have on each line

##################
# Advanced - combustion flamelet parameters (Recommend using default values unless you know what you are doing)
##################
# if you wanted to use custom values set the variable Combustion_Flamelet_Default to False
Combustion_Flamelet_Default = False
Combustion_Flamelet_Fourier_Number_Initial = 1
Combustion_Flamelet_Fourier_Number_Increase_Factor = 1.5
Combustion_Flamelet_ODE_Relative_Error_Tolerance = 1e-05
Combustion_Flamelet_ODE_Absolute_Error_Tolerance = 1e-15
Combustion_Flamelet_Convergence_Tolerance = 1e-05
Combustion_Flamelet_Maximum_Integration_Time = 1000
Combustion_Flamelet_Scalar_Dissipation_Initial = 0.01
Combustion_Flamelet_Scalar_Dissipation_Multiplier = 10
Combustion_Flamelet_Scalar_Dissipation_Step = 5
Combustion_Flamelet_Grid_Points_Initial = 8
Combustion_Flamelet_Grid_Points_Max = 64
Combustion_Flamelet_Value_Ratio_Max_Change = 0.5
Combustion_Flamelet_Slope_Ratio_Max_Change = 0.5
Combustion_Flamelet_Max_number = 10

##################
# Advanced - combustion PDF parameters (Recommend using default values unless you know what you are doing)
##################
# if you wanted to use custom values set the variable Combustion_PDF_Default to False
Combustion_PDF_Default = False
Combustion_PDF_Grid_Points_Initial = 15 # Default = 15
Combustion_PDF_Grid_Points_Max = 200 # Default = 200
Combustion_PDF_Value_Ratio_Max_Change = 0.25 # Default = 0.25
Combustion_PDF_Slope_Ratio_Max_Change = 0.25 # Default = 0.25
Combustion_PDF_Species_Max = 20 # Default = 20
Combustion_PDF_Tempurature_Min = 200 # Default = 298, in Kelvin
Combustion_PDF_Max_Heat_Loss_Fraction = 100 # Default = 0.667
Combustion_PDF_Max_Heat_Gain_Fraction = 1 # Default = 0.25


################################################################################################################################################################################
# Some input checks and auto fixes enabled
################################################################################################################################################################################
# Auto disable GPU solver if combustion sim is nessary
if Solver_Run_Combustion == True:
    Solver_GPU_Enhanced = False

if Combustion_Flamelet_Default == True:
    # Flamelet parameters are set to default therefore parameters are being corrected to default values
    Combustion_Flamelet_Fourier_Number_Initial = 1
    Combustion_Flamelet_Fourier_Number_Increase_Factor = 2
    Combustion_Flamelet_ODE_Relative_Error_Tolerance = 1e-05
    Combustion_Flamelet_ODE_Absolute_Error_Tolerance = 1e-15
    Combustion_Flamelet_Convergence_Tolerance = 1e-05
    Combustion_Flamelet_Maximum_Integration_Time = 1000
    Combustion_Flamelet_Scalar_Dissipation_Initial = 0.01
    Combustion_Flamelet_Scalar_Dissipation_Multiplier = 10
    Combustion_Flamelet_Scalar_Dissipation_Step = 5
    Combustion_Flamelet_Grid_Points_Initial = 8
    Combustion_Flamelet_Grid_Points_Max = 64
    Combustion_Flamelet_Value_Ratio_Max_Change = 0.5
    Combustion_Flamelet_Slope_Ratio_Max_Change = 0.5
    Combustion_Flamelet_Max_number = 10

if Combustion_PDF_Default == True:
    # PDF parameters are set to default therefore parameters are being corrected to default values
    Combustion_PDF_Grid_Points_Initial = 15
    Combustion_PDF_Grid_Points_Max = 200
    Combustion_PDF_Value_Ratio_Max_Change = 0.25
    Combustion_PDF_Slope_Ratio_Max_Change = 0.25
    Combustion_PDF_Species_Max = 20
    Combustion_PDF_Tempurature_Min = 298
    Combustion_PDF_Max_Heat_Loss_Fraction = 0.667
    Combustion_PDF_Max_Heat_Gain_Fraction = 0.25

if Combustion_Fuel_Species_Name == "ch4":
    # Setting variable for the carbon number of fuel for nox simulation
    Combustion_Fuel_Carbon_Number = 1
else:
    # Fuel is nh3 or h2 so no carbon
    Combustion_Fuel_Carbon_Number = 0

# Adding reports of species
if Combustion_Report_Fuel == True:
    # Adding fuel species to report
    Combustion_Report_Species += [Combustion_Fuel_Species_Name]

if Combustion_Report_Carbons == True:
    # Adding c, co and co2 to report
    Combustion_Report_Species += ["c","co","co2"]

if Combustion_Report_NOx == True:
    if Combustion_Fuel_Species_Name == "nh3":
        Combustion_Report_Species += ["no","no2","nh2"]
    else:
        Combustion_Report_Species += ["no"]

# Adding contour plots of species
if Combustion_Contour_Plot_Fuel == True:
    # Adding fuel species to contor plots
    Combustion_Contour_Species_To_Plot += [Combustion_Fuel_Species_Name]

if Combustion_Contour_Plot_NOx == True:
    if Combustion_Fuel_Species_Name == "nh3":
        Combustion_Contour_Species_To_Plot += ["no","no2","nh2"]
    else:
        Combustion_Contour_Species_To_Plot += ["no"]

if Combustion_Contour_Plot_Carbons == True:
    # Adding contour plots of carbons
    Combustion_Contour_Species_To_Plot += ["c","co","co2"]

################################################################################################################################################################################
# Function Definitions
################################################################################################################################################################################
# Generic Linear interpolator calculator
# For a known value of x, a value of y can be found using two known points of x and y, one before and after x.
def Function_Linear_Interpolator(x0,y0,x1,y1,x):
    y = y0 + (((x-x0)*(y1-y0))/(x1-x0))
    return(y)

# Generic Turbulent Intensity Calculator
# For use in the turbulant physics models
def Function_Turbulent_Intensity(
                                Mass_Flow_Rate,
                                Fluid_Density,
                                Charactoristic_Length,
                                Area,
                                Dynamic_Viscosity
                                ):
    
    Fluid_Velocity = Mass_Flow_Rate / (Area * Fluid_Density)
    Re = (Fluid_Density * Fluid_Velocity * Charactoristic_Length)/(Dynamic_Viscosity)

    return (0.16 * pow(Re,-0.125))

# Function to pull the dimensions out of the file name string formate bellow
# a<Value For a>b<Value For b>l<Value for l>IDAM<Value for IDAM>IDAI<Value for IDAI>AMNUM<Value for AMNUM>
def Function_Get_Dimensions_From_String(Input_File_Name):

    a = float(Input_File_Name[Input_File_Name.find("a") + 1:Input_File_Name.find("b")])
    b = float(Input_File_Name[Input_File_Name.find("b") + 1:Input_File_Name.find("l")])
    l = float(Input_File_Name[Input_File_Name.find("l") + 1:Input_File_Name.find("IDAM")])
    idam = float(Input_File_Name[Input_File_Name.find("IDAM") + 4:Input_File_Name.find("IDAI")])
    idai = float(Input_File_Name[Input_File_Name.find("IDAI") + 4:Input_File_Name.find("AMNUM")])
    amnum = float(Input_File_Name[Input_File_Name.find("AMNUM") + 5:Input_File_Name.find("FC")])
    fc = float(Input_File_Name[Input_File_Name.find("FC") + 2:Input_File_Name.find("EL")])
    el = float(Input_File_Name[Input_File_Name.find("EL") + 2:])

    return(a,b,l,idam,idai,amnum,fc,el)

# Function to find the mesh sizes based on the table lookup
def Function_Cylinder_Mesh_Size_Finder(Table_Mesh_Size_Tube_Radius, Table_Mesh_Size, Given_Radius):
    # First check to see if tube radius given is already in table.
    if Table_Mesh_Size_Tube_Radius.count(Given_Radius) > 0: 
        Provided_Mesh_Size = Table_Mesh_Size[Table_Mesh_Size_Tube_Radius.index(Given_Radius)]
    else:
        Temp_Index_Below_Given = 0
        Temp_Index_Above_Given = len(Table_Mesh_Size_Tube_Radius) - 1
        # Radius is inbetween values in the radius mesh table
        for i in range(len(Table_Mesh_Size_Tube_Radius)):
            if Table_Mesh_Size_Tube_Radius[i] < Given_Radius and Table_Mesh_Size_Tube_Radius[i] > Table_Mesh_Size_Tube_Radius[Temp_Index_Below_Given]:
                Temp_Index_Below_Given = i

            elif  Table_Mesh_Size_Tube_Radius[i] > Given_Radius and Table_Mesh_Size_Tube_Radius[i] < Table_Mesh_Size_Tube_Radius[Temp_Index_Above_Given]:
                Temp_Index_Above_Given = i
        
        # Now find the mesh size using linear interpolation
        Provided_Mesh_Size = Function_Linear_Interpolator(
                                                        Table_Mesh_Size_Tube_Radius[Temp_Index_Below_Given],
                                                        Table_Mesh_Size[Temp_Index_Below_Given],
                                                        Table_Mesh_Size_Tube_Radius[Temp_Index_Above_Given],
                                                        Table_Mesh_Size[Temp_Index_Above_Given],
                                                        Given_Radius
                                                        )

    return (Provided_Mesh_Size)

# Function to get x cordinate for monitoring points
def Function_Get_X_Val(Radius, Number_Of_Tubes):
    Angle = math.radians((360 / (Number_Of_Tubes * 2)))
    
    x = math.cos(Angle) * Radius
    
    return (x)

# Function to get y cordinate for monitoring points
def Function_Get_Y_Val(Radius, Number_Of_Tubes):
    Angle = math.radians((360 / (Number_Of_Tubes * 2)))

    y = -1 * math.sin(Angle) * Radius

    return (y)


################################################################################################################################################################################
# Initializing for loop and finding common values to use in meshing and solution
################################################################################################################################################################################
for i in List_File_Names_To_Run:
    Parameter_A, Parameter_B, Parameter_L, Parameter_IDAM, Parameter_IDAI, Parameter_AMNUM, Parameter_FC, Parameter_EL = Function_Get_Dimensions_From_String(i)
    
    # default setup for meshing and solver parameters
    Meshing_Acceptable = False
    Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides = Default_Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides


    ################################################################################################################################################################################
    # Automated Meshing Section
    ################################################################################################################################################################################
    while Meshing_Run == True and Meshing_Acceptable == False:

        # Launching Meshing workflow
        Meshing_Session = pyfluent.launch_fluent(
                                                version="3d",
                                                precision="double",
                                                processor_count=4,
                                                mode="meshing",
                                                show_gui= Meshing_Show_GUI,
                                                cwd= (Parent_Folder + "/" + i)
                                                )
        Meshing_Session.health_check_service.is_serving

        # Setup for workflow and meshing variables/ classes to simplfy later code
        workflow = Meshing_Session.workflow
        meshing = Meshing_Session.meshing

        #######################################################
        # Initializing Workflow and setting units
        #######################################################
        workflow.InitializeWorkflow(WorkflowType='Watertight Geometry')
        meshing.GlobalSettings.LengthUnit.set_state('mm')
        meshing.GlobalSettings.AreaUnit.set_state('mm^2')
        meshing.GlobalSettings.VolumeUnit.set_state('mm^3')

        #######################################################
        # Importing the model to the mesher
        #######################################################
        workflow.TaskObject['Import Geometry'].Arguments.set_state({
                                                                    'FileName': (Parent_Folder + "/" + i + "/" + i + ".dsco"),
                                                                    'ImportCadPreferences': {'ShowImportCadPreferences': False}
                                                                })
        workflow.TaskObject['Import Geometry'].Execute()

        #######################################################
        # Adding Local Sizes to Model (Mesh refinement)
        #######################################################
        # Air Inlet Refinements
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Air_Inlet_Side_Refinement',
                                                'BOIFaceLabelList': ['air_inlet_side'],
                                                'BOISize': Function_Cylinder_Mesh_Size_Finder(
                                                                                            Table_Of_standard_Mesh_Sizes_Tube_Radius,
                                                                                            Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides,
                                                                                            (Parameter_IDAI/2)
                                                                                            )
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Ammonia Inlet Refinements
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Ammonia_Inlet_Side_Refinement',
                                                'BOIFaceLabelList': ['ammonia_inlet_sides'],
                                                'BOISize': Function_Cylinder_Mesh_Size_Finder(
                                                                                            Table_Of_standard_Mesh_Sizes_Tube_Radius,
                                                                                            Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides,
                                                                                            (Parameter_IDAM/2)
                                                                                            )
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Exaust Refinements
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Outlet_Side_Refinement',
                                                'BOIFaceLabelList': ['outlet_side'],
                                                'BOISize': Function_Cylinder_Mesh_Size_Finder(
                                                                                            Table_Of_standard_Mesh_Sizes_Tube_Radius,
                                                                                            Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides,
                                                                                            (Parameter_B)
                                                                                            )
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Combustor Refinement
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Combustor_Side_Refinement',
                                                'BOIFaceLabelList': ['combustor_side'],
                                                'BOISize': Function_Cylinder_Mesh_Size_Finder(
                                                                                            Table_Of_standard_Mesh_Sizes_Tube_Radius,
                                                                                            Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides,
                                                                                            (Parameter_A)
                                                                                            )
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Periodic Boundry faces refinement
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Periodic_Boundry_Faces_Refinement',
                                                'BOIFaceLabelList': ['periodic_boundry_sidea','periodic_boundry_sideb'],
                                                'BOISize': Meshing_Periodic_Boundry_Element_Size
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Combustor Fillets Refinement
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Combustor_Fillet',
                                                'BOIFaceLabelList': ['combustor_fillet'],
                                                'BOISize': Function_Cylinder_Mesh_Size_Finder(
                                                                                            Table_Of_standard_Mesh_Sizes_Tube_Radius,
                                                                                            Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides,
                                                                                            (Parameter_FC)
                                                                                            )
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        # Outlet face refinement
        workflow.TaskObject['Add Local Sizing'].Arguments.set_state({
                                                'AddChild': 'yes',
                                                'BOIControlName': 'Combustor_Outlet_Face_Refinement',
                                                'BOIFaceLabelList': ['outlet_face'],
                                                'BOISize': Meshing_Combustor_Outlet_Face_Element_Size
                                                            })
        workflow.TaskObject['Add Local Sizing'].AddChildAndUpdate()

        #######################################################
        # Generate the Surface Mesh
        #######################################################
        workflow.TaskObject['Generate the Surface Mesh'].Arguments.set_state({
                                                                            'CFDSurfaceMeshControls': {'CellsPerGap': 2,'CurvatureNormalAngle': 10},
                                                                            'SurfaceMeshPreferences': {'ShowSurfaceMeshPreferences': True}
                                                                            })
        workflow.TaskObject['Generate the Surface Mesh'].Execute()

        #######################################################
        # Add in the Periodic Boundry Condition
        #######################################################
        workflow.TaskObject['Generate the Surface Mesh'].InsertNextTask(CommandName='SetUpPeriodicBoundaries')
        workflow.TaskObject['Set Up Periodic Boundaries'].Arguments.set_state({
                                                                                'LCSOrigin': {'OriginX': 0,'OriginY': 0,'OriginZ': 0},
                                                                                'LCSVector': {'VectorX': 0,'VectorY': 0,'VectorZ': 1},
                                                                                'LabelList': ['periodic_boundry_sidea', 'periodic_boundry_sideb'],
                                                                                'ListAllLabelToggle': False,
                                                                                'MeshObject': '',
                                                                                'Method': 'Automatic - pick both sides',
                                                                                'PeriodicityAngle': (360/Parameter_AMNUM),
                                                                                'RemeshBoundariesOption': 'yes',
                                                                                'SelectionType': 'label',
                                                                                'TransShift': {'ShiftX': 0,'ShiftY': 0,'ShiftZ': 1},
                                                                                'Type': 'Rotational'
                                                                            })
        workflow.TaskObject['Set Up Periodic Boundaries'].Execute()

        #######################################################
        # Describe the Geometry
        #######################################################
        workflow.TaskObject['Describe Geometry'].UpdateChildTasks(SetupTypeChanged=False)
        workflow.TaskObject['Describe Geometry'].Arguments.set_state({
                                                                    'SetupType': 'The geometry consists of only fluid regions with no voids'
                                                                    })
        workflow.TaskObject['Describe Geometry'].UpdateChildTasks(SetupTypeChanged=True)
        workflow.TaskObject['Describe Geometry'].Execute()

        #######################################################
        # Update Boundries and Regions
        #######################################################
        workflow.TaskObject['Update Boundaries'].Arguments.set_state({
                                                                    'BoundaryLabelList': [
                                                                                        'air_inlet_face',
                                                                                        'air_inlet_side',
                                                                                        'ammonia_inlet_sides',
                                                                                        'ammonia_inlet_faces',
                                                                                        'outlet_side'
                                                                                        ],
                                                                    'BoundaryLabelTypeList': [
                                                                                            'mass-flow-inlet',
                                                                                            'wall',
                                                                                            'wall',
                                                                                            'mass-flow-inlet',
                                                                                            'wall'
                                                                                            ],
                                                                    'OldBoundaryLabelList': [
                                                                                            'air_inlet_face',
                                                                                            'air_inlet_side',
                                                                                            'ammonia_inlet_sides',
                                                                                            'ammonia_inlet_faces',
                                                                                            'outlet_side'
                                                                                            ],
                                                                    'OldBoundaryLabelTypeList': [
                                                                                                'velocity-inlet',
                                                                                                'velocity-inlet',
                                                                                                'velocity-inlet',
                                                                                                'velocity-inlet',
                                                                                                'pressure-outlet'
                                                                                                ],
                                                                    })
        workflow.TaskObject['Update Boundaries'].Execute()
        workflow.TaskObject['Update Regions'].Execute()

        #######################################################
        # Adding Boundry Layers
        #######################################################
        workflow.TaskObject['Add Boundary Layers'].Arguments.set_state({
                                                                        'BLControlName': 'Sides_Boundry_Layer',
                                                                        'LocalPrismPreferences': {'ShowLocalPrismPreferences': False},
                                                                        'NumberOfLayers': Boundry_Layer_Settings_Number_Of_layers,
                                                                        'FirstHeight': Boundry_Layer_Settings_First_Height,
                                                                        'OffsetMethodType': 'uniform',
                                                                        'Rate': Boundry_Layer_Settings_Rate,
                                                                        'ZoneSelectionList': [
                                                                                            'air_inlet_side',
                                                                                            'ammonia_inlet_sides',
                                                                                            'combustor_side',
                                                                                            'combustor_faces',
                                                                                            'combustor_fillet',
                                                                                            'outlet_side'
                                                                                            ]
                                                                    })
        workflow.TaskObject['Add Boundary Layers'].AddChildAndUpdate()

        #######################################################
        # Generate the Volme Mesh
        #######################################################
        workflow.TaskObject['Generate the Volume Mesh'].Arguments.set_state({
                                                                            'PrismPreferences': {r'ShowPrismPreferences': False},
                                                                            'VolumeMeshPreferences': {'ShowVolumeMeshPreferences': False}
                                                                            })
        workflow.TaskObject['Generate the Volume Mesh'].Execute()

        #######################################################
        # Perform Mesh Check, Saving File and ending session or relooping untill element count is bellow 1 million
        #######################################################
        Meshing_Session.tui.mesh.check_mesh()
        Meshing_Session.tui.mesh.check_quality()
        print("Add refiement inputs here make some adjustments here")
        Temp_Input = input("Is the mesh acceptable? (Type y (lower case) for yes, Type refine (lower case) to decrease element size by: " + str(Meshing_Element_Size_Decrease_Ammount) + "mm, all other values will remesh by increasing element size by: " + str(Meshing_Element_Size_Increase_Ammount) + ". : ")

        if Temp_Input == "y":
            # Mesh has been termed as acceptable and the mesh will be saved.
            Meshing_Acceptable = True
            Meshing_Session.tui.file.write_case(i + "_Mesh")
            # Note that the above saves to the current workspace location

        elif Temp_Input == "refine":
            # Mesh is not acceptable (element count too high) remeshing with larger element sizes
            for iterable in range(len(Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides)):
                Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides[iterable] -= Meshing_Element_Size_Decrease_Ammount
      

        else:
            # Mesh is not acceptable (element count too high) remeshing with larger element sizes
            for iterable in range(len(Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides)):
                Table_Of_standard_Mesh_Sizes_Element_Size_Cylinder_Sides[iterable] += Meshing_Element_Size_Increase_Ammount
      
        # Closing Session either for new geometry or for remeshing
        Meshing_Session.exit()


    ################################################################################################################################################################################
    # Automated Solver Section
    ################################################################################################################################################################################
    if Solver_Run == True:
        if Solver_GPU_Enhanced == True:
            # Launching FLuent Solver Session with gpu acceleration (Need Nividia GPU see Ansys Documentation for compatable GPUs)
            Solver_Session = pyfluent.launch_fluent(
                                                    version="3d",
                                                    precision="double",
                                                    processor_count=3,
                                                    gpu=True,
                                                    mode="solver",
                                                    show_gui= Solver_Show_GUI,
                                                    cwd= (Parent_Folder + "/" + i)
                                                )
        
        else:
            # Launching Fluent without GPU acceleration only using CPU compute
            Solver_Session = pyfluent.launch_fluent(
                                                    version="3d",
                                                    precision="double",
                                                    processor_count=4,
                                                    mode="solver",
                                                    show_gui= Solver_Show_GUI,
                                                    cwd= (Parent_Folder + "/" + i)
                                                )
            
        #######################################################
        # Importing Mesh to solver and performing mesh check
        #######################################################
        Solver_Session.file.read(file_type="case", file_name=(Parent_Folder + "/" + i + "/" + i + "_Mesh" + ".cas.h5"))
        Solver_Session.tui.mesh.check()
        Solver_Session.tui.define.units("length", "mm")

        #######################################################
        # Model parameters setup depending on turbulance only or combustion and turbulance
        #######################################################
        # Changing the turbulance model to ke realizable
        Solver_Session.tui.define.models.viscous.ke_realizable("yes")
        
        if Solver_Run_Combustion == True:
            ##############################################################################################################
            # Turbulant Combustion specific setup
            ##############################################################################################################

            # Energy Model Enable
            Solver_Session.setup.models.energy.enabled = True

            # Non Premixed Combustion TUI command and answers
            # Note that some are editable in this script but some are not to make this script work.
            Solver_Session.tui.define.models.species.non_premixed_combustion(
                                                                            "yes", # Enable non premix combustion model?
                                                                            "yes", # Create the PDF table?
                                                                            "no", # Adiabatic?
                                                                            "no", # Equilibrium (else Laminar Flamelet)?
                                                                            "yes", # Steady Diffusion Flamelet (else Unsteady)?
                                                                            "yes", # Create Flamelet (else Import)?
                                                                            "yes", # Automatic Refinement?
                                                                            Combustion_Operating_Pressure, # Operating Pressure (in [Pa])
                                                                            ('"' + Parent_Folder + "/" + Combustion_Kinetics_Thermodynamics_Tranport_Folder_Name +"/" + Combustion_File_Name_Kinetics + '"'), # CHEMKIN Mechanism File location path
                                                                            ('"' + Parent_Folder + "/" + Combustion_Kinetics_Thermodynamics_Tranport_Folder_Name +"/" + Combustion_File_Name_Thermodynamics + '"'), # Thermodynamic Data File location path
                                                                            Fuel_Inlet_Tempurature, # Fuel Stream Temperature
                                                                            Air_Inlet_Tempurature, # Oxidizer Stream Temperature (Air in this combustion sim)
                                                                            "yes", # Specify Species in mass fractions (next set of numbers based on mass)
                                                                            ('"' + Combustion_Fuel_Species_Name + '"'), # Fuel species name (Note that you could do mixtures of different fuels but, functionality wold have to be added into the script)
                                                                            1, # Fuel Fraction of specified species
                                                                            0, # Oxidizer Fraction of specified species
                                                                            ('"' + "o2" + '"'), # Species Name
                                                                            0, # Fuel Fraction of specified species
                                                                            0.233, # Oxidizer Fraction of specified species (specifically oxygen component of air)
                                                                            ('"' + "n2" + '"'), # Species Name
                                                                            0, # Fuel Fraction of specified species
                                                                            0.767, # Oxidizer Fraction of specified species (specifically nitrogen component of air)
                                                                            '""', # To confirm all species of interest have been added to simulation
                                                                            Combustion_Flamelet_Fourier_Number_Initial, # Initial Fourier Number
                                                                            Combustion_Flamelet_Fourier_Number_Increase_Factor, # Fourier Number Increase Factor
                                                                            Combustion_Flamelet_ODE_Relative_Error_Tolerance, # ODE Relative Error Tolerance
                                                                            Combustion_Flamelet_ODE_Absolute_Error_Tolerance, # ODE Absolute Error Tolerance
                                                                            Combustion_Flamelet_Convergence_Tolerance, # Flamelet Convergence Tolerance
                                                                            Combustion_Flamelet_Maximum_Integration_Time, # Maximum Integration Time
                                                                            Combustion_Flamelet_Scalar_Dissipation_Initial, # Initial Scalar Dissipation
                                                                            Combustion_Flamelet_Scalar_Dissipation_Multiplier, # Scalar Dissipation Multiplier
                                                                            Combustion_Flamelet_Scalar_Dissipation_Step, # Scalar Dissipation Step
                                                                            Combustion_Flamelet_Grid_Points_Initial, # Initial Number of Grid Points in Flamelet
                                                                            Combustion_Flamelet_Grid_Points_Max, # Maximum Number of Grid Points in Flamelet
                                                                            Combustion_Flamelet_Value_Ratio_Max_Change, # Maximum Change in Value Ratio
                                                                            Combustion_Flamelet_Slope_Ratio_Max_Change, # Maximum Change in Slope Ratio
                                                                            Combustion_Flamelet_Max_number, # Maximum Number of Flamelets
                                                                            "no", # User defined flamelet parameters?
                                                                            "yes", # Write Flamelet file?
                                                                            ('"' + Parent_Folder + "/" + i + "/" + i + "_Flamelet.fla" + '"'), # Flamelet File Path to save to
                                                                            "yes", # Automatic PDF Grid Refinement?
                                                                            Combustion_PDF_Grid_Points_Initial, # Initial Number of Grid Points
                                                                            Combustion_PDF_Grid_Points_Max, # Maximum Number of Grid Points
                                                                            Combustion_PDF_Value_Ratio_Max_Change, # Maximum Change in Value Ratio
                                                                            Combustion_PDF_Slope_Ratio_Max_Change, # Maximum Change in Slope Ratio
                                                                            Combustion_PDF_Species_Max, # Maximum Number of Species
                                                                            Combustion_PDF_Tempurature_Min, # Minimum Temperature
                                                                            "yes", # Include Equilibrium Flamelet?
                                                                            "no", # Ok to discard PDF file?
                                                                            ('"' + Parent_Folder + "/" + i + "/" + i + "_PDF.pdf" + '"'), # PDF file path to save to
                                                                            "no" # Write in binary format?
                                                                            )

            # Changing Non Premixed expert parameters for better PDF table generation
            Solver_Session.tui.define.models.species.non_premixed_combustion_expert(
                                                                                1, # Probability Density Function for PDF Table (0 for double-delta 1 for beta)
                                                                                "yes", # Second order interpolation for PDF Table (else 4th order)
                                                                                1, # Maximum Temperature to Adiabatic Flame Temperature Ratio
                                                                                Combustion_PDF_Max_Heat_Loss_Fraction, # Maximum Heat Loss Fraction
                                                                                Combustion_PDF_Max_Heat_Gain_Fraction, # Maximum Heat Gain Fraction
                                                                                "yes", # Enable checking of PDF table temperature limits?
                                                                                1, # Equilibrium Calculations at Inlets: 0:full-equilibrium, 1:frozen-composition
                                                                                1 # Verbosity
                                                                                )

            # Enabling the new expert parameters for the pdf table
            Solver_Session.tui.define.models.species.non_premixed_combustion(
                                                                            "yes", # Enable non premix combustion model?
                                                                            "yes", # Create the PDF table?
                                                                            "no", # Adiabatic?
                                                                            "no", # Equilibrium (else Laminar Flamelet)?
                                                                            "yes", # Steady Diffusion Flamelet (else Unsteady)?
                                                                            "no", # Create Flamelet (else Import)?
                                                                            "yes", # Standard Flamelet?
                                                                            ('"' + Parent_Folder + "/" + i + "/" + i + "_Flamelet.fla" + '"'), # Flamelet File Path to import
                                                                            '""', # Exiting menu to choose flamets to import
                                                                            "yes", # Automatic PDF Grid Refinement?
                                                                            Combustion_PDF_Grid_Points_Initial, # Initial Number of Grid Points
                                                                            Combustion_PDF_Grid_Points_Max, # Maximum Number of Grid Points
                                                                            Combustion_PDF_Value_Ratio_Max_Change, # Maximum Change in Value Ratio
                                                                            Combustion_PDF_Slope_Ratio_Max_Change, # Maximum Change in Slope Ratio
                                                                            Combustion_PDF_Species_Max, # Maximum Number of Species
                                                                            Combustion_PDF_Tempurature_Min, # Minimum Temperature
                                                                            "yes", # Include Equilibrium Flamelet?
                                                                            "no", # Ok to discard PDF file?
                                                                            ('"' + Parent_Folder + "/" + i + "/" + i + "_PDF.pdf" + '"'), # PDF file path to save to
                                                                            "no" # Write in binary format?
                                                                            )

            #######################################################
            # NOx Parameter adjustments
            #######################################################
            if Combustion_NOx_Enable == True:
                # Enabling NOx modeling in combustion
                Solver_Session.tui.define.models.nox("yes")

                # Enabling NOx chemistry and changing some parameters
                Solver_Session.tui.define.models.nox_parameters.nox_chemistry(
                                                                            "yes", # Thermal NOx?
                                                                            Combustion_NOx_O_Rad, # Model Choice for O radicals (1 oxygen atom): 0 = equilibrium, 1 = partial-equilibrium, 2 = instantaneous
                                                                            Combustion_NOx_OH_Rad, # Model Choice for OH radicals (1 oxygen and hydrogen atom): 0 = equilibrium, 1 = partial-equilibrium, 2 = instantaneous
                                                                            "yes", # Prompt NOx?
                                                                            ('"' + Combustion_Fuel_Species_Name + '"'), # Fuel Species Name
                                                                            '""', # To exit fuel species stream
                                                                            "no", # Remove Fuel Species from list?
                                                                            Combustion_Fuel_Carbon_Number, # Fuel carbon number
                                                                            Combustion_Equivalance_Ratio, # Equivalence ratio
                                                                            "no", # N2O Path?
                                                                            "no", # Include reburn?
                                                                            "no"  # Include SNCR?
                                                                            )

                # Adjusting NOx parameters
                Solver_Session.tui.define.models.nox_parameters.nox_turbulence_interaction(
                                                                                        "yes", # Enable turbulance interation
                                                                                        1, # PDF mode options: 1 = tempurature, 4 = mixture fraction. Script only setup to deal with tempurature mode
                                                                                        Combustion_PDF_Species_Max, # PDF points
                                                                                        0, # PDF Type mode option: 0 = beta, 1 = gaussian. Script only setup for beta
                                                                                        1, # Tempurature variance options: 0 = algebraic, 1 = transported. Script only setup for transported option
                                                                                        0, # Tmax options: 0 = global tmax, 1 = local tmax factor, 2 = specified tmax. Script  setup for global tmax option
                                                                                        "yes" # Use cpu intensive (accurate) option for pdf integration?
                                                                                        )
        
            #######################################################
            # Boundry conditions
            #######################################################
            # Air Inlet boundry condition setup
            Air_Inlet = Solver_Session.setup.boundary_conditions.mass_flow_inlet["air_inlet_face"]
            Air_Inlet.mass_flow = Air_Total_Flow_Rate / Parameter_AMNUM
            Air_Inlet.ke_spec = "Intensity and Hydraulic Diameter"
            Air_Inlet.turb_intensity =  Function_Turbulent_Intensity(
                                                                    Air_Total_Flow_Rate,
                                                                    Air_Density,
                                                                    Parameter_IDAI,
                                                                    (math.pi * ((Parameter_IDAI/2) ** 2)),
                                                                    Air_Dynamic_Viscosity
                                                                    )
            Air_Inlet.turb_hydraulic_diam = float(Parameter_IDAI)
            Air_Inlet.t0 = Air_Inlet_Tempurature
            Air_Inlet.fmean = 0 # Mixture fraction, as just air set to 0 as 1 is pure fuel

            # Ammonia Inlet boundry condition setup
            Ammonia_Inlet = Solver_Session.setup.boundary_conditions.mass_flow_inlet["ammonia_inlet_faces"]
            Ammonia_Inlet.mass_flow = Fuel_Total_Flow_Rate / Parameter_AMNUM
            Ammonia_Inlet.ke_spec = "Intensity and Hydraulic Diameter"
            Ammonia_Inlet.turb_intensity =  Function_Turbulent_Intensity(
                                                                        (Fuel_Total_Flow_Rate / Parameter_AMNUM),
                                                                        Fuel_Denisty,
                                                                        Parameter_IDAM,
                                                                        (math.pi * ((Parameter_IDAM/2) ** 2)),
                                                                        Fuel_Dynamic_Viscosity
                                                                        )
            Ammonia_Inlet.turb_hydraulic_diam = float(Parameter_IDAM)
            Ammonia_Inlet.t0 = Fuel_Inlet_Tempurature
            Ammonia_Inlet.fmean = 1 # Mixture fraction, as just fuel set to 1 as 1 is pure fuel

            # Pressure Outlet
            Outlet_Face = Solver_Session.setup.boundary_conditions.pressure_outlet["outlet_face"]

            #######################################################
            # Setting up tempurature convergence criteria
            #######################################################
            Solver_Session.tui.solve.report_definitions.add(
                                                            "Tempurature_Report", # Report definition name
                                                            "surface-massavg", # Type of report
                                                            "field", # Opening Field Sepection promt
                                                            "temperature", # Selecting Tempurature Field
                                                            "surface-names", # Selecting the Surface
                                                            "outlet_face", # The Surface Name to use
                                                            "" # Exit report definition
                                                            )
            

            # Report file creation
            Solver_Session.tui.solve.report_files.add(
                                                    "Tempurature_Report", # Report file Name
                                                    "report-defs", # Linking report definition to the file
                                                    "Tempurature_Report", # Report definition to link
                                                    "" # To exit report definition linking
                                                    )
                
            # Report plot creation
            Solver_Session.tui.solve.report_plots.add(
                                                    "Tempurature_Report", # Report plot Name
                                                    "report-defs", # Linking report definition to the plot
                                                    "Tempurature_Report", # Report definition to link
                                                    "" # To exit report definition linking
                                                    )


            Solver_Session.tui.solve.convergence_conditions(
                                                            "conv-reports", # Going into convergence reports options
                                                            "add", # Add a convergence criteria
                                                            "Outlet_Tempurature_Average", # Name of convergence report
                                                            "initial-values-to-ignore", # Number of values to ignore at first
                                                            20, # Number of values to ignore at first
                                                            "previous-values-to-consider", # Number of previous values to see how much they change by
                                                            15, # Number of previous values to see how much they change by
                                                            "print?", # Print to the terminal?
                                                            "yes", # Print to the erminal yes option selected
                                                            "report-defs", # Selecting the related report file to look at
                                                            "Tempurature_Report", # Report Name
                                                            "stop-criterion", # Defining the stop criterion
                                                            Combustion_Convergence_Tempurature, # Defining the difference value
                                                            )

            #######################################################
            # Setting up extra reports of species at the output face
            #######################################################
            for Species_To_Report in Combustion_Report_Species:
                # Report definitions setup
                Solver_Session.tui.solve.report_definitions.add(
                                                                (Species_To_Report + "_Out"), # Report Definition Name
                                                                "surface-massavg", # Type of report
                                                                "field", # Enable choosing field
                                                                (Species_To_Report), # Species option
                                                                "surface-names", # Choosing surface name that this will be conputed from
                                                                "outlet_face", # selected surface to compute from
                                                                "" # To close surface selection
                                                                )
                
                # Report file creation
                Solver_Session.tui.solve.report_files.add(
                                                        (Species_To_Report + "_Out"), # Report file Name
                                                        "report-defs", # Linking report definition to the file
                                                        (Species_To_Report + "_Out"), # Report definition to link
                                                        "" # To exit report definition linking
                                                        )
                
                # Report plot creation
                Solver_Session.tui.solve.report_plots.add(
                                                        (Species_To_Report + "_Out"), # Report plot Name
                                                        "report-defs", # Linking report definition to the plot
                                                        (Species_To_Report + "_Out"), # Report definition to link
                                                        "" # To exit report definition linking
                                                        )

                # Convergence Criteria For respective species added
                Solver_Session.tui.solve.convergence_conditions(
                                                                "conv-reports", # Going into convergence reports options
                                                                "add", # Add a convergence criteria
                                                                (Species_To_Report + "_Converge"), # Name of convergence report
                                                                "initial-values-to-ignore", # Number of values to ignore at first
                                                                20, # Number of values to ignore at first
                                                                "previous-values-to-consider", # Number of previous values to see how much they change by
                                                                15, # Number of previous values to see how much they change by
                                                                "print?", # Print to the terminal?
                                                                "yes", # Print to the terminal yes option selected
                                                                "report-defs", # Selecting the related report file to look at
                                                                (Species_To_Report + "_Out"), # Report Name
                                                                "stop-criterion", # Defining the stop criterion
                                                                Combustion_Report_Refinement_Value # Defining the difference value
                                                                )

            #######################################################
            # Setting up monitoring points in the combustor
            #######################################################
            if Combustion_Monitor_Points_Enable == True:
                for Points in range(Combustion_Monitor_Points_Ammount_In_Line): # Range defines how many points will be in in the combustor length
                    # Adding points in line with outlet tube radius
                    Length_Along_Line_Temp = ((Parameter_L / (Combustion_Monitor_Points_Ammount_In_Line + 1)) * (Points + 1))
                    Solver_Session.tui.surface.point_surface(
                                                            ("Point_Outlet_" + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Point Name Given
                                                            Function_Get_X_Val(Parameter_B, Parameter_AMNUM), # x (in [mm]) (Math done to get point on radius of outlet)
                                                            Function_Get_Y_Val(Parameter_B, Parameter_AMNUM), # y (in [mm]) (Math done to get point on radius of outlet)
                                                            Length_Along_Line_Temp # z (in [mm]) (Is the dimensiosn of combustor length, Note, base of combustor is on the xy plane)
                                                            )

                    # Adding points in line with inlet tube radius
                    Solver_Session.tui.surface.point_surface(
                                                            ("Point_Inlet_" + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Point Name Given
                                                            Function_Get_X_Val((Parameter_IDAI / 2), Parameter_AMNUM), # x (in [mm]) (Math done to get point on radius of inlet)
                                                            Function_Get_Y_Val((Parameter_IDAI / 2), Parameter_AMNUM), # y (in [mm]) (Math done to get point on radius of inlet)
                                                            Length_Along_Line_Temp # z (in [mm]) (Is the dimensiosn of combustor length, Note, base of combustor is on the xy plane)
                                                            )
                    
                    # Adding points in line with half inlet tube radius
                    Solver_Session.tui.surface.point_surface(
                                                            ("Point_Half_Inlet_" + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Point Name Given
                                                            Function_Get_X_Val((Parameter_IDAI / 4), Parameter_AMNUM), # x (in [mm]) (Math done to get point on radius of half inlet)
                                                            Function_Get_Y_Val((Parameter_IDAI / 4), Parameter_AMNUM), # y (in [mm]) (Math done to get point on radius of half inlet)
                                                            Length_Along_Line_Temp # z (in [mm]) (Is the dimensiosn of combustor length, Note, base of combustor is on the xy plane)
                                                            )

                    # Adding monitors and convergence criteria for all points
                    for Report in ["Outlet_","Inlet_","Half_Inlet_"]:
                        #######################################################
                        # Adding monitors and convergence criteria for tempurature of all points
                        #######################################################
                        # Generating Report definition
                        Solver_Session.tui.solve.report_definitions.add(
                                                                        ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition name
                                                                        "surface-vertexavg", # Type of report
                                                                        "field", # Opening Field Sepection promt
                                                                        "temperature", # Selecting Tempurature Field
                                                                        "surface-names", # Selecting the Surface
                                                                        ("Point_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # The Surface Name to use
                                                                        "" # Exit report definition
                                                                        )
                        
                        # Report file creation
                        Solver_Session.tui.solve.report_files.add(
                                                                ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report file Name
                                                                "report-defs", # Linking report definition to the file
                                                                ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition to link
                                                                "" # To exit report definition linking
                                                                )
                            
                        # Report plot creation
                        Solver_Session.tui.solve.report_plots.add(
                                                                ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report plot Name
                                                                "report-defs", # Linking report definition to the plot
                                                                ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition to link
                                                                "" # To exit report definition linking
                                                                )

                        # Adding in convergence criteria
                        Solver_Session.tui.solve.convergence_conditions(
                                                                        "conv-reports", # Going into convergence reports options
                                                                        "add", # Add a convergence criteria
                                                                        ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Name of convergence report
                                                                        "initial-values-to-ignore", # Number of values to ignore at first
                                                                        20, # Number of values to ignore at first
                                                                        "previous-values-to-consider", # Number of previous values to see how much they change by
                                                                        15, # Number of previous values to see how much they change by
                                                                        "print?", # Print to the terminal?
                                                                        "yes", # Print to the erminal yes option selected
                                                                        "report-defs", # Selecting the related report file to look at
                                                                        ("Tempurature_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report Name
                                                                        "stop-criterion", # Defining the stop criterion
                                                                        Combustion_Monitor_Points_Convergence_Number # Defining the difference value
                                                                        )
                        
                        #######################################################
                        # Adding monitors and convergence criteria for velocity headed towards outlet at all points
                        #######################################################
                        # Generating Report definition
                        Solver_Session.tui.solve.report_definitions.add(
                                                                        ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition name
                                                                        "surface-vertexavg", # Type of report
                                                                        "field", # Opening Field Sepection promt
                                                                        "z-velocity", # Selecting Tempurature Field
                                                                        "surface-names", # Selecting the Surface
                                                                        ("Point_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # The Surface Name to use
                                                                        "" # Exit report definition
                                                                        )
                        
                        # Report file creation
                        Solver_Session.tui.solve.report_files.add(
                                                                ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report file Name
                                                                "report-defs", # Linking report definition to the file
                                                                ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition to link
                                                                "" # To exit report definition linking
                                                                )
                            
                        # Report plot creation
                        Solver_Session.tui.solve.report_plots.add(
                                                                ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report plot Name
                                                                "report-defs", # Linking report definition to the plot
                                                                ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report definition to link
                                                                "" # To exit report definition linking
                                                                )

                        # Adding in convergence criteria
                        Solver_Session.tui.solve.convergence_conditions(
                                                                        "conv-reports", # Going into convergence reports options
                                                                        "add", # Add a convergence criteria
                                                                        ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Name of convergence report
                                                                        "initial-values-to-ignore", # Number of values to ignore at first
                                                                        20, # Number of values to ignore at first
                                                                        "previous-values-to-consider", # Number of previous values to see how much they change by
                                                                        15, # Number of previous values to see how much they change by
                                                                        "print?", # Print to the terminal?
                                                                        "yes", # Print to the erminal yes option selected
                                                                        "report-defs", # Selecting the related report file to look at
                                                                        ("Flow_Report_" + Report + str(int(Length_Along_Line_Temp)) + "_Combustor_Base"), # Report Name
                                                                        "stop-criterion", # Defining the stop criterion
                                                                        Combustion_Monitor_Points_Convergence_Number # Defining the difference value
                                                                        )


        else:
            ##############################################################################################################
            # Turbulance model only specific setup
            ##############################################################################################################
            # Add the relevent Materials to model (ie ammonia)
            Solver_Session.setup.materials.database.copy_by_name(type="fluid", name="ammonia-vapor")

            #######################################################
            # Boundry conditions
            #######################################################
            # Air Inlet boundry condition setup
            Air_Inlet = Solver_Session.setup.boundary_conditions.mass_flow_inlet["air_inlet_face"]
            Air_Inlet.mass_flow = Air_Total_Flow_Rate / Parameter_AMNUM
            Air_Inlet.ke_spec = "Intensity and Hydraulic Diameter"
            Air_Inlet.turb_intensity =  Function_Turbulent_Intensity(
                                                                    Air_Total_Flow_Rate,
                                                                    Air_Density,
                                                                    Parameter_IDAI,
                                                                    (math.pi * ((Parameter_IDAI/2) ** 2)),
                                                                    Air_Dynamic_Viscosity
                                                                    )
            Air_Inlet.turb_hydraulic_diam = Parameter_IDAI

            # Ammonia Inlet boundry condition setup
            Ammonia_Inlet = Solver_Session.setup.boundary_conditions.mass_flow_inlet["ammonia_inlet_faces"]
            Ammonia_Inlet.mass_flow = Fuel_Total_Flow_Rate / Parameter_AMNUM
            Ammonia_Inlet.ke_spec = "Intensity and Hydraulic Diameter"
            Ammonia_Inlet.turb_intensity =  Function_Turbulent_Intensity(
                                                                        (Fuel_Total_Flow_Rate / Parameter_AMNUM),
                                                                        Fuel_Denisty,
                                                                        Parameter_IDAM,
                                                                        (math.pi * ((Parameter_IDAM/2) ** 2)),
                                                                        Fuel_Dynamic_Viscosity
                                                                        )
            Ammonia_Inlet.turb_hydraulic_diam = Parameter_IDAM

            # Pressure Outlet
            Outlet_Face = Solver_Session.setup.boundary_conditions.pressure_outlet["outlet_face"]

        #######################################################
        # Coupled Pressure solver options
        #######################################################

        # Adjusting settings for coupled pressure solver
        if Solver_GPU_Enhanced == False and Solver_Run_Combustion == True:

            # Enabling the coupled pressure velcity solver.
            Solver_Session.tui.solve.set.p_v_coupling(24) # Options are as follows from Ansys Manual: 20 = Simple, 21 = SimpleC, 22 = PISO, 24 = Coupled, 25 = Fractional Step


            Solver_Session.tui.solve.set.pseudo_time_method.global_time_step_settings(
                                                                                    "yes",
                                                                                    1,
                                                                                    Solver_Session_Coupled_Solver_Time_Scale_Factor_Global
                                                                                    )

            # For Energy 
            Solver_Session.tui.solve.set.pseudo_time_method.advanced_options(
                                                                            "enthalpy",
                                                                            "yes",
                                                                            Solver_Session_Coupled_Solver_Time_Scale_Factor_Energy
                                                                            )

            # For f mean fraction
            Solver_Session.tui.solve.set.pseudo_time_method.advanced_options(
                                                                            "fmean",
                                                                            "yes",
                                                                            Solver_Session_Coupled_Solver_Time_Scale_Factor_Mean_Mixture_Fraction
                                                                            )

        #######################################################
        # Initializing and running solver
        #######################################################
        # Disabeling residual plotting
        #Solver_Session.tui.solve.monitors.residual.plot("no")

        Solver_Session.solution.initialization.hybrid_initialize()

        Solver_Session.solution.run_calculation.iterate(iter_count=Solver_Max_Iteration_Count)

        #######################################################
        # Post Processing Plots for turbulant and combustion simulations
        #######################################################
        # Generating mesh for mesh visualisation
        Solver_Session.tui.display.objects.create(
                                                "mesh", # type of graphic to create
                                                "Combustor_Mesh", # Name of new graphic
                                                "surface-list", # Opening surface selection
                                                "air_inlet_face", # Selecting Surface
                                                "air_inlet_side", # Selecting Surface
                                                "combustor_faces", # Selecting Surface
                                                "combustor_side", # Selecting Surface
                                                "periodic_boundry_sidea", # Selecting Surface
                                                "periodic_boundry_sideb", # Selecting Surface
                                                "ammonia_inlet_faces", # Selecting Surface
                                                "ammonia_inlet_sides", # Selecting Surface
                                                "outlet_side", # Selecting Surface
                                                "outlet_face", # Selecting Surface
                                                "" # Closing Surface selection
                                                )

        # Generating pathlines from inlet faces described above to export in simulation report
        for Pathlines_To_Display in Solver_Pathlines_Faces_Output:
            Solver_Session.tui.display.objects.create(
                                                    "pathlines", # Type of graphic to create
                                                    ('"Pathlines_From_' + Pathlines_To_Display + '"'), # Name of new graphic
                                                    "step", # Ammount of steps pathlines will take
                                                    Solver_Pathlines_Steps, # Ammount of steps
                                                    "surfaces-list", # Opening surface selection
                                                    Pathlines_To_Display, # Selecting Surface
                                                    "()", # Closing Surface selection
                                                    "field", # Entering field selection menu
                                                    "velocity-magnitude" # selecting field to display
                                                    )

        if Solver_Run_Combustion == True:
            #######################################################
            # Post Processing Plots for combustion simulations only
            #######################################################
            if Combustion_Contour_Plot_Tempurature == True:
                # Making contor of tempurature
                Solver_Session.tui.display.objects.create(
                                                        "contour", # Type of graphic to create
                                                        "Tempurature_Contor_Plot", # Name of contour graphic
                                                        "contour-lines?", # Asking to enable contour lines
                                                        "yes", # Enabling contour lines
                                                        "field", # Opening menu to choose field to display
                                                        "temperature", # Choosing field to display
                                                        "surfaces-list", # Opening menu to select surface to dsiplay on
                                                        "periodic_boundry_sidea", # Selected surface
                                                        "()" # Closing surface selection
                                                        )

            if Combustion_Contour_Plot_Species == True:
                for Contour_Species_To_Display in Combustion_Contour_Species_To_Plot:
                    Solver_Session.tui.display.objects.create(
                                                            "contour", # Type of graphic to create
                                                            (Contour_Species_To_Display + "_Contour_Plot"), # Contour plot species
                                                            "contour-lines?", # Asking to enable contour lines
                                                            "yes", # Enabling contour lines
                                                            "field", # Opening menu to choose field to display
                                                            Contour_Species_To_Display, # Choosing spcies to display
                                                            "surfaces-list", # Opening menu to select surface to dsiplay on
                                                            "periodic_boundry_sidea", # Selected surface
                                                            "()" # Closing surface selection
                                                            )

        #######################################################
        # Generating and exporting simulation report to HTML and  pdf form
        #######################################################
        # Generating simulation report file
        Solver_Session.tui.report.simulation_reports.generate_simulation_report('"' + i + "_Sim_Data" + '"')

        # exporting simulation report to pdf and HTML formats
        Solver_Session.tui.report.simulation_reports.export_simulation_report_as_html(
                                                                                    ('"' + i + "_Sim_Data" + '"'), # Report name to export
                                                                                    ('"' + Parent_Folder + "/" + i + "/" + i + '_Sim_Data_HTML"')# File path to save simulation report
                                                                                    )
        
        Solver_Session.tui.report.simulation_reports.export_simulation_report_as_pdf(
                                                                                    ('"' + i + "_Sim_Data" + '"'), # Report name to export
                                                                                    ('"' + Parent_Folder + "/" + i + "/" + i + '_Sim_Data_PDF"')# File path to save simulation report
                                                                                    )

        #######################################################
        # Saving case and data file and closing session
        #######################################################
        Solver_Session.file.write(file_type = "case", file_name = (i + "_Result.cas.h5"))
        Solver_Session.file.write(file_type = "data", file_name = (i + "_Result.dat.h5"))

        Solver_Session.exit()

        time.sleep(60) # sleep point to allow fluent session to close properly for next fluent session to load correctly.

