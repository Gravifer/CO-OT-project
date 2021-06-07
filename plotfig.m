% iter = 7000;
% sxx = 1:1:iter;
% plot(sxx(1:iter),(valds(1:iter)),'b');
% %hold on
% %plot(sxx,valds,'b');
% %legend('the prime objective value','the dual objective value');
% xlabel('number of iteration','FontSize',35); 
% ylabel('the dual objective value','FontSize',35,'rotation',90)
% %ylabel('the prime objective value','FontSize',35,'rotation',90)
% %ylabel('log(error)','FontSize',35,'rotation',90)
% title('The convergence of admmdual','FontSize',35)


sxx = 1:1:iter;
plot(sxx(1:iter), sinkhorn_opts(1:iter),'b');
% hold on
% plot(sxx(1:iter),sinkhorn_opts(1:iter),'b');
% legend('ipot','sinkhorn');
xlabel('number of iteration','FontSize',35); 
ylabel('the objective value','FontSize',35,'rotation',90)
% title('The convergence of sinkhorn','FontSize',35)
title('The convergence of sinkhorn','FontSize',35)