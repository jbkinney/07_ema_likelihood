clear all
close all

% Collect muk results
directory = 'abf1_muk/'
meanEmats = [];
for k=[1,2,3,4,5,6,8,9,10]
    file_name = [directory 'run_' num2str(k)];
    disp(file_name)
    load(file_name)
    N = numel(r.models);
    emats = zeros(r.L,4,N);
    for n=1:N
        model = fixModelGauge(r.models(n));
        emats(:,:,n) = model.emat;
    end
    emats(:,:,1:1000) = []; % Remove first 1000 models
    meanEmats(:,:,k) = mean(emats,3);
end
figure
meanEmats(:,:,7) = [];
for i=1:9
    subplot(5,2,i)
    imagesc(meanEmats(:,:,i)')
end
abf1_muk_emat(:,:) = mean(meanEmats,3);


% Collect lee results
directory = 'abf1_lee/'
meanEmats = [];
for k=1:10
    file_name = [directory 'run_' num2str(k)];
    disp(file_name)
    load(file_name)
    N = numel(r.models);
    emats = zeros(r.L,4,N);
    for n=1:N
        model = fixModelGauge(r.models(n));
        emats(:,:,n) = model.emat;
    end
    emats(:,:,1:1000) = []; % Remove first 1000 models
    meanEmats(:,:,k) = mean(emats,3);
end
figure
for i=1:10
    subplot(5,2,i)
    imagesc(meanEmats(:,:,i)')
end
abf1_lee_emat(:,:) = mean(meanEmats,3);

save abf1_mean_emats.mat abf1_muk_emat abf1_lee_emat