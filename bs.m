% band-stop filter
% input : data(ch,dp)->both of them will work, band1->49, band2->51 etc? 
% Output : filtered data

function vargout = bs(data,fs,jisu,band1,band2)

% band-stop
[b,a] = butter(jisu,[(band1)/(fs/2),(band2)/(fs/2)],'stop');
if length(data(:,1)) > length(data(1,:))
    vargout = filtfilt(b,a,data);
else
    data2 = data';
    data2 = filtfilt(b,a,data2);
    vargout = data2';
end

end