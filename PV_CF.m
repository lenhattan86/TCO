Location = {'12-1-2011','12-2-2011','12-3-2011','12-4-2011','12-5-2011','12-6-2011','12-7-2011'};
for i = 1:1:size(Location,2)
    pv = load(['C:\Users\liuzhenh.AMERICAS\Downloads\traces\PV\papvdata-',char(Location(i)),'.txt']);
    if size(pv,1) < 280
        'Error'
    end
    pvm(i) = mean(pv(:,2))*12/130;
end
pvm


Location = {'3-1-2011','3-2-2011','3-3-2011','3-4-2011','3-5-2011','3-6-2011','3-7-2011'};
for i = 1:1:size(Location,2)
    pv = load(['C:\Users\liuzhenh.AMERICAS\Downloads\traces\PV\papvdata-',char(Location(i)),'.txt']);
    if size(pv,1) < 280
        'Error'
    end
    pvm(i) = mean(pv(:,2))*12/130;
end 
pvm


Location = {'7-1-2011','7-2-2011','7-3-2011','7-4-2011','7-5-2011','7-6-2011','7-7-2011'};

for i = 1:1:size(Location,2)
    pv = load(['C:\Users\liuzhenh.AMERICAS\Downloads\traces\PV\papvdata-',char(Location(i)),'.txt']);
    pvm(i) = mean(pv(:,2))*12/130;
end 


Location = {'9-1-2011','9-2-2011','9-3-2011','9-4-2011','9-5-2011','9-6-2011','9-7-2011'};
for i = 1:1:size(Location,2)
    pv = load(['C:\Users\liuzhenh.AMERICAS\Downloads\traces\PV\papvdata-',char(Location(i)),'.txt']);
    if size(pv,1) < 280
        'Error'
    end
    pvm(i) = mean(pv(:,2))*12/130;
end 
pvm


