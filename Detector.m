classdef Detector < handle
 %Detector осуществляет обнаружение целей

 properties
 %определяем параметры детектора
 thr %порог обнаружения
 Fa %вероятность ложных тревог

 detections %выходной массив данных
 end

 methods
 function obj = Detector(fa)
 obj.Fa = fa;
 end
 function detect(obj,inData,r,a) %обнаружение целей
 %вычисляем порог в заданной области свободной от целей r - a
 temp = abs(inData(r,a));
 s2 = var(temp(:));
 obj.thr = 2.5*sqrt(-2*log(obj.Fa)*s2);
 obj.detections = inData > obj.thr;
 end
 function show(obj,hInData)
 figure;
 %сырые данные
 subplot(1,2,1);
 temp = abs(hInData(:));
 q = quantile(temp,0.9);
 imagesc(abs(hInData),[0 q]);
 xlabel('Угол, отс');
 ylabel('Дальность, отс');

 if ~isempty(obj.detections)
 subplot(1,2,2);
 %данные после порога
 imagesc(abs(obj.detections),[0 1]);
 xlabel('Угол, отс');
 ylabel('Дальность, отс');
 end
 end
 end
end