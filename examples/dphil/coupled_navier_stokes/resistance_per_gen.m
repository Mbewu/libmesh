file_input = importdata('results/resistance_per_generation.dat');
resistance = file_input.data;

for i=1:length(resistance(:,1))
    resistance(i,3) =  resistance(i,2)/2.^(i-1);
end

figure
plot(resistance(:,1),resistance(:,2))
%plot(resistance(:,1),resistance(:,3))