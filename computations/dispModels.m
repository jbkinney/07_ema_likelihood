function dispModels(models)
    imageFileLocation = 'Characters/';
    bases = ['ACGT'];   
    baseImages = {};
    for i=1:length(bases)
        fname = [imageFileLocation bases(i) '.png'];
        baseImages{i} = imread(fname);
    end
    
    cla % clear axes
    M = size(models,1);
    N = size(models,2);
    set(gca, 'YLim', M*[0 2.01], 'XLim', N*[0 1])
    numTrials = 1e6/(M*N);
    for m=1:M
        for n=1:N
            % Displays the sequence logo of an energy matrix
            model = models(m,n);
            pwm = model2pwm(model, numTrials);
            L = size(pwm,1);
            heights = 2 + sum((pwm.*log2(pwm + eps))')';
            
            
            for i=1:L
                charheight = 0;
                B= sortrows([pwm(i,:)' (1:4)'], [1]);
                baseindices =  B(:,2);
                charheights = heights(i)*B(:,1);
                chartops = cumsum(charheights')';
                charbots = chartops-charheights;
                for j=1:4
                    basepic = baseImages{baseindices(j)};
                    top = chartops(j);
                    bot = charbots(j);
                    if top-bot > .001
                        image('CData', basepic,...
                        'XData', (n-1)+[i-1 i]/L, 'YData', 2*(m-1)+[top bot]); 
                    end
                end
            end
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
        end
    end
    for m=0:M
        line([0 N], 2*[m m], 'Color', 'k', 'linewidth', 1);
    end
    for n=0:N
        line([n n], [0 2*M], 'Color', 'k', 'linewidth', 1);
    end
    set(gca, 'XColor', 'w')
    set(gca, 'YColor', 'w')
    line(N*[0 1], 2*M*[0 0], 'Color', 'k', 'linewidth', 1);
    line(N*[0 1], 2*M*[1 1], 'Color', 'k', 'linewidth', 1);
    line(N*[0 0], 2*M*[0 1], 'Color', 'k', 'linewidth', 1);
    line(N*[1 1], 2*M*[0 1], 'Color', 'k', 'linewidth', 1);      
end