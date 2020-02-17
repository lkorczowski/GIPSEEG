function freqs=freqs_val(Fs,nb_points)
%freqs=freqs_val(Fs,nb_points)

Vec=1:nb_points;
nygf=(Fs);
freqs=0:nygf/nb_points:nygf-nygf/nb_points;