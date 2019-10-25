function data = resample_phins(data)

att = zeros(3,size(data,1));
for i=1:size(data,1)
R = [data(i,20),data(i,21),data(i,22);data(i,23),data(i,24),data(i,25);data(i,26),data(i,27),data(i,28)];
att(:,i) = rot2rph(R);
end

att = resample2(data(1:5000/10:end-5000/10,12),att(:,5000/10:5000/10:end)',data(:,12),'linear');

for i=1:size(data,1)

    R = rph2R(att(i,:));
    data(i,20:28) = reshape(R',1,9);
    
end
