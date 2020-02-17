#!/usr/bin/env ruby
require "openvibe-launcher"

DESIGNER_PATH = 'C:/openvibe/openvibe-designer.cmd'
SCENARIO_PATH = 'F:/BciTmp/P300-Marco/Data-Old-Marco/csv_converter.xml'
MODE = "--play-fast" 
OPTIONS = ["--no-gui","--no-session-management"]

Openvibe = Openvibe::Launcher.new(DESIGNER_PATH,SCENARIO_PATH,MODE,OPTIONS)

searchpath = File.join("**", "*.ov")

files = Dir.glob(searchpath).each do |path|
  env = Hash.new
  destfile = path.gsub(".ov",".csv")
  env["sourcefile"] = File.absolute_path(path)
  env["destfile"] = File.absolute_path(destfile)
  puts env["sourcefile"]
  #if !(File.exist?(destfile))
    Openvibe.setEnv(env)
	Openvibe.setOutput("NUL")
	Openvibe.start
  #end
end
