function textBox(str)
% function textBox(str)
% put small box right bottom for mandatory information in plot
dim = [.01 .01 .02 .02];
%str = ['bi2013a s' num2str(tix) ' C'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','linestyle','none');
end