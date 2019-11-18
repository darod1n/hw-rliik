classdef AKP < handle
 %AKP осуществляет автокомпенсацию активных помех

 properties
 %определяем параметры автокомпенсации
 wnd = 30 %число отсчетов для оценки коэфф. компенсации
 g = 10 %коэффициент компенсации

 th_shift = 5000; %смещение для заднего лепестка 5000 = -180 гр
 k_back = 10^(-13.3/20); %коэффициент усиления компенсационной антенны

 akp_data %выходной массив данных
 akp_ang %выходной массив кодов углов данных
 end

 methods
 function compensate(obj,hInData,hInAng) %компенсация данных
 obj.akp_data = []; %удаляем результаты прошлых вычислений
 obj.akp_ang = [];

 %вычисляем сигнал действующий на входе
 a = (hInAng - obj.th_shift)>1;
 [~,col] = find(a > 0,1,'first');
 for ii = col:length(hInAng)
 u0 = hInData(:,ii); %Формула 20
 angIdx = ii-obj.th_shift:ii-1;
 angIdx(angIdx < 1) = [];
 u1 = sum(hInData(:,angIdx),2)*obj.k_back; %формула 21 и 24
 K = obj.g*sum(u0(end-obj.wnd:end).*u1(end-obj.wnd:end),1)/...
 (1+obj.g*sum(u1(end-obj.wnd:end).^2,1));%формула 23
 obj.akp_data = [obj.akp_data (u0 - K*u1)];%Формула 22
 obj.akp_ang = [obj.akp_ang hInAng(ii)];
 end
 end

 function show(obj,hInData,hInAng)
 figure;
 %сырые данные
 subplot(1,2,1);
 temp = abs(hInData(:));
 q = quantile(temp,0.9);
 idxAng = (hInAng - obj.th_shift) > 1;
 imagesc(abs(hInData(:,idxAng)),[0 q]);
 xlabel('Угол, отс');
 ylabel('Дальность, отс');

 subplot(1,2,2);
 %данные после компенсации
 imagesc(abs(obj.akp_data),[0 q]);
 xlabel('Угол, отс');
 ylabel('Дальность, отс');
 end

 function Kp = calc_kp(obj,r,a,hInData)
 %оценка коэффициента подавления в выбранной точке
 Kp = 20*log10(abs(obj.akp_out)/abs(hInData(r,a)));
 end
 end
end