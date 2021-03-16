clearvars
close all

klist = linspace(-pi,pi,200);
t = 0;

eig_val = zeros(2,length(klist));
eig_vec1 = zeros(2,length(klist));
eig_vec2 = zeros(2,length(klist));

for i = 1:length(klist)
    k = klist(i);
    H = [[-2*cos(k),-2*cos(k)-t];[-2*cos(k)-t,-2*cos(k)]];
    [V,D] = eig(H);
    
    eig_val(:,i) = diag(D);
    eig_vec1(:,i) = V(:,1);
    eig_vec2(:,i) = V(:,2);
end

hold on
plot(klist/pi,eig_val(1,:))
plot(klist/pi,eig_val(2,:))
hold off

