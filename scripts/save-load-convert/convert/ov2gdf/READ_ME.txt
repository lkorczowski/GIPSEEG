# This ruby script will convert all the openvibe files (.ov) into .gdf files
# it use the openvibe scenario gdf_converter.xml
# please modify this scenario if you like with openvibe-designer (by default, there is a temporal filer)
#
#
#0. you need to have openvibe installed
#1. you need to install ruby from https://www.ruby-lang.org/fr/ 
#2. you need to install openvibe-launcher gem from https://rubygems.org/gems/openvibe-launcher
#3. then run the ruby script with command ruby for exemple:
# run cmd.exe and put in the command line >D:\convert_all.rb
#
# please read each comment inside convert_all.rb  in order to configure the parameters
# when the script has finished properly, you'll have a message. DO NOT CLOSE the cmd or the data could be corrupted.
# if you have closed the script before it is over, kill all the openvibe related processus in the task manager such as:
# cmd.exe, conhost.exe, openvibe-designer.exe, etc.
#
# IMPORTANT :
# if you want to run openvibe in background, you need to open openvibe-designer.cmd and modify 
#SET OV_PAUSE=PAUSE
#to
#SET OV_PAUSE=
#without this setting, you need to press a button every time the script has converted a file
#
# ***********Last modification**********************
# Louis Korczowski, march 2015 @GIPSA-lab
# **************************************************
