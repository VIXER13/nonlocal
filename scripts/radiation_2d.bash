#Use this command before explotation:
#sed -i 's/\r$//' radiation_2d.bash

line="------------------------------------------------------"
echo "Hello! You are working with thermal conductivity equation solver."
echo $line
echo "Searching for build directory..."
build_dir=../build/
if [ -d $build_dir ]; then
	echo "The $build_dir directory exists."
	cd $build_dir
else
	echo "The $build_dir directory does not exist"
	ls
fi
echo $line
echo "Starting make thermal_nonstationary_radiation_2d... "
time make thermal_nonstationary_radiation_2d 
echo $line
main_dir=build/calculations
way_to_folder="D:/nonlocal/"$main_dir
function calculation {
file_name="$6"
echo $file_name 
mkdir $way_to_folder/$2/$file_name  
# Input format [program name] <path to mesh> <r1> <r2> <p1> <save_path>
time ./src/2d/thermal/thermal_nonstationary_radiation_2d $1 $3 $4 $5 $way_to_folder/$2/$file_name 
}
echo "Starting calculations..."
echo "All results will be written in folder $way_to_folder ."
echo $line

#------------------------------------------------------------------------------------------------------------------------
# КАК ПРОВОДИТЬ ЕДИНИЧНЫЙ РАСЧЕТ:
# 1) создание папки для сохранения данных : mkdir $way_to_folder/<название_папки>
# 2) calculation <путь к сетке> <папка расчета>  <r1> <r2> <p1> <название расчета>
#------------------------------------------------------------------------------------------------------------------------


mkdir $way_to_folder/$folder
calculation  "D:/nonlocal/meshs/mesh_rectangular_2d.su2" "tests" 1 1 1 "test2d_nonlocal"
