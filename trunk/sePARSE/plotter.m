%% plotter
% a script to hold all of the messy plot/drawing stuff and get it out of
% the main code.
% use: whichplot='psis';plotter;fprintf('%s',whichplot)
if (plot_a.show) 
    set(0,'CurrentFigure',plot_a.fig); % instead of figure(plot_a) to avoid focus grabbing
    switch whichplot
        case 'gradients' % show parameter gradients                   
            subplot(3,3,1);
            plot(h(Voxels),'.');            % M grad
            hold on
            plot(h(Voxels+numVox),'.r');    % exponent grad
            plot(h(Voxels+numVox*2),'.g');  % R prime
            hold off
        case 'slice' % show an example "slice" of the trajectory
            subplot(3,3,7);
            plot(kss(1:N1));
        case 'progress' % display the images of each parameter to show progress
            set(gcf,'Name',[num2str(iteration) ' iterations']);
            subplot(3,3,8);
            semilogy((CC(1:iteration))); 
            title([num2str(timeLength) '  --  ' num2str(CC(iteration))]);
            subplot(3,3,9);
            plot(frRect(:,10:10:50))
            subplot(3,3,2);
            imagesc(frRect);
            axis image;
            title('f');
            subplot(3,3,5);
            imagesc(M0Rect);
            axis image;
            title('M0');
            subplot(3,3,3);
            imagesc(R2Rect);
            axis image;
            title('R2');
            subplot(3,3,6);
            imagesc(R2primeRect);
            axis image;
            title('R2prime');
            colormap gray;
        otherwise
    end
    drawnow
    whichplot=''; 
else
    whichplot='.';
end
