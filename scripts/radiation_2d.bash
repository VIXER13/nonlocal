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
time ./src/2d/thermal/thermal_nonstationary_radiation_2d $1 $3 $4 $5 $way_to_folder/$2/$file_name $7
}
echo "Starting calculations..."
echo "All results will be written in folder $way_to_folder ."
echo $line

#------------------------------------------------------------------------------------------------------------------------
# КАК ПРОВОДИТЬ ЕДИНИЧНЫЙ РАСЧЕТ:
# 1) создание папки для сохранения данных : mkdir $way_to_folder/<название_папки>
# 2) calculation <путь к сетке> <папка расчета>  <r1> <r2> <p1> <название расчета>
#------------------------------------------------------------------------------------------------------------------------


mkdir $way_to_folder/"calcs" 
calculation  "D:/nonlocal/meshs/25x25.su2" "calcs" 0.8 0.8 0.5 "25x25_nonlocal_0_01" 0.01
calculation  "D:/nonlocal/meshs/50x50.su2" "calcs" 0.8 0.8 0.5 "50x50_nonlocal_0_01" 0.01
calculation  "D:/nonlocal/meshs/5x5.su2" "calcs" 0.8 0.8 0.5 "5x5_nonlocal_0_01" 0.01
calculation  "D:/nonlocal/meshs/0_02.su2" "calcs" 0.8 0.8 0.5 "0_02_nonlocal_0_01" 0.01
calculation  "D:/nonlocal/meshs/0_04.su2" "calcs" 0.8 0.8 0.5 "0_04_nonlocal_0_01" 0.01

calculation  "D:/nonlocal/meshs/25x25.su2" "calcs" 0.8 0.8 0.5 "25x25_nonlocal_0_005" 0.005
calculation  "D:/nonlocal/meshs/50x50.su2" "calcs" 0.8 0.8 0.5 "50x50_nonlocal_0_005"  0.005
calculation  "D:/nonlocal/meshs/5x5.su2" "calcs" 0.8 0.8 0.5 "5x5_nonlocal_0_005"  0.005
calculation  "D:/nonlocal/meshs/0_02.su2" "calcs" 0.8 0.8 0.5 "0_02_nonlocal_0_005" 0.005
calculation  "D:/nonlocal/meshs/0_04.su2" "calcs" 0.8 0.8 0.5 "0_04_nonlocal_0_005"  0.005

calculation  "D:/nonlocal/meshs/25x25.su2" "calcs" 0.8 0.8 0.5 "25x25_nonlocal_0_0025" 0.0025
calculation  "D:/nonlocal/meshs/50x50.su2" "calcs" 0.8 0.8 0.5 "50x50_nonlocal_0_0025" 0.0025
calculation  "D:/nonlocal/meshs/5x5.su2" "calcs" 0.8 0.8 0.5 "5x5_nonlocal_0_0025" 0.0025
calculation  "D:/nonlocal/meshs/0_02.su2" "calcs" 0.8 0.8 0.5 "0_02_nonlocal_0_0025" 0.0025
calculation  "D:/nonlocal/meshs/0_04.su2" "calcs" 0.8 0.8 0.5 "0_04_nonlocal_0_0025" 0.0025
