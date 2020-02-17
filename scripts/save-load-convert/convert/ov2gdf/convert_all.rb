#!/usr/bin/env ruby
# D:\data\convert_all.rb
#
#1. you need to install ruby from https://www.ruby-lang.org/fr/ 
#2. install https://rubygems.org/pages/download
#2. you need to install openvibe-launcher gem from https://rubygems.org/gems/openvibe-launcher
#3. then launch this file with command ruby for exemple:
#open cmd.exe and put in the command line >D:\convert_all.rb
#
# please read each comment inside this code in order to configure the parameters
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

require "openvibe-launcher"
puts "hello"
DESIGNER_PATH = 'D:/openvibe/dist/openvibe-designer.cmd'
#SCENARIO_PATH = 'D:/Mes Documents GIPSA/MATLAB/LK_TOOLBOX/scripts/save-load-convert/convert/ov2gdf/gdf_converter.xml'
SCENARIO_PATH = 'D:/data/erwanRAW/gdf_converter.xml'
#DATA_PATH='D:/Mes Documents GIPSA/MATLAB/LK_TOOLBOX/data/MARTI/' #put the path where the script will search and convert ALL the .ov files (can be long)
DATA_PATH='D:/data/Hyperscanning/MARTI/EEG Screening BI march-april 2015/'
MODE = "--play-fast" 
OPTIONS = ["--no-gui","--no-session-management"]
#OPTIONS = ["--no-session-management"]
#OPTIONS = ["--run-bg"]

Openvibe = Openvibe::Launcher.new(DESIGNER_PATH,SCENARIO_PATH,MODE,OPTIONS)
searchpath = File.join(DATA_PATH,'**', "*.ov")
compt=0;
files = Dir.glob(searchpath).each do |path|
env = Hash.new
destfile = path.gsub(".ov",".gdf")
env["sourcefile"] = File.absolute_path(path)
env["destfile"] = File.absolute_path(destfile)
puts env["sourcefile"]
#if !(File.exist?(env["destfile"])) #uncomment this line if you don't want to override the existing .gdf files
puts env["destfile"]
Openvibe.setEnv(env)
Openvibe.setOutput("a") #"NUL" for window, "/dev/null" for linux
Openvibe.start
puts "file written"
compt=compt+1;
#if compt==2
#puts "function aborded premarturly by the user"
#puts "Program finished, #{compt} files written"
#return
#end
#else
#puts "file already exists"
#end
  end
	puts "Program finished, #{compt} files written"
