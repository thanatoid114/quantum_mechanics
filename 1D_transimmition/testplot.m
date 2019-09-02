clear all
close all

load parameters
disp('making plot');
NT=length(time);

v = VideoWriter('mov.avi');
open(v);


for k = 1:100:NT
    plot(x,abs(psi(k,:).^2),'b');
    ylim([0 1.2]);
    %     plot(x,abs(fftshift(fft(fftshift(psi(k,:))))),'b-o');
    title(['Time: ' sprintf('%1.5f',time(1,k)) ' Sum: ' sprintf('%2.5f',sum(abs(psi(k,:)))) ]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);
disp('finish');