#!/usr/bin/env ruby
require "openvibe-launcher"

DESIGNER_PATH = '/home/kirsh/openvibe/trunk/dist/openvibe-designer.sh'
SCENARIO_PATH = '/home/kirsh/Documents/Data/Erwan/Data/gdf_converter.xml'
MODE = "--play-fast" 
OPTIONS = ["--no-gui","--no-session-management"]

Openvibe = Openvibe::Launcher.new(DESIGNER_PATH,SCENARIO_PATH,MODE,OPTIONS)

searchpath = File.join("**", "*.ov")

files = Dir.glob(searchpath).each do |path|
  env = Hash.new
  destfile = path.gsub(".ov",".gdf")
  env["sourcefile"] = File.absolute_path(path)
  env["destfile"] = File.absolute_path(destfile)
  puts env["sourcefile"]
  if !(File.exist?(destfile))
    Openvibe.setEnv(env)
    Openvibe.setOutput("/dev/null")
    Openvibe.start
  end
end
