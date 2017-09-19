function draw_MPc(varargin)
% 	
% draw_MPc
% draw_MPc(renorm): renormalization
% draw_MPc(renorm,colors):  colors must be a 4x3 matrix containing the
% color vectors for the metal, carbon, nitrogen and hydrogen
% 

nargs = nargin;

%% get axis data

% oldax=gca;
% newax=axes('position',get(gca,'position'),'color','none','visible','off','xlim',get(oldax,'xlim'),'ylim',get(oldax,'ylim'),'zlim',get(oldax,'zlim'));
% newax=axes('position',get(gca,'position'),'color','none','visible','off','xlim',get(oldax,'xlim'),'ylim',get(oldax,'ylim'),'zlim',[-1 1]);

%% Farben der Atome: y m c r g b w k

FarbeC_0	= [.2 .2 .2];        %[.5 .5 .5]
FarbeH_0	= 'b';
FarbeN_0	= 'g';
FarbeCu_0	= 'y';

if nargs ==0
	renorm  = 1;	
	FarbeC  = FarbeC_0;
	FarbeH  = FarbeH_0;
	FarbeN  = FarbeN_0;
	FarbeCu = FarbeCu_0;
elseif nargs == 1;
	renorm = varargin{1};
	FarbeC  = FarbeC_0;
	FarbeH  = FarbeH_0;
	FarbeN  = FarbeN_0;
	FarbeCu = FarbeCu_0;
elseif nargs >=2;
	renorm = varargin{1};
	FarbeC  = 'c';varargin{2}(2,:);
	FarbeH  = varargin{2}(4,:);
	FarbeN  = varargin{2}(3,:);
	FarbeCu = varargin{2}(1,:);
end

C = coordinates('Pc');

%% Größen der Atome; Cu=135 pm, N=65pm, C=70pm, H=25pm 

scaling = .5;

rC  = scaling*0.7;
rH  = scaling*0.25;
rN  = scaling*0.65;
rCu = scaling*1.;
Color = 'k';
Linewidth = 1.5;



hold on


%% Molekülgitter dazumalen:


A = C(1:14,:);         % unit cell Koordinaten
R = [0 -1;1 0];             % Drehmatrix

for k=1:4
    B = (R*A')';                        % für Verbindung zwischen Atom 14 und Atom 2 der nächsten Unit Cell
    plot( A(1:4,1) , A(1:4,2) ,'Color',Color,'LineWidth',Linewidth)         % die ersten 4 Verbindungen


    for i=0:3
        plot( [A(4+2*i , 1) A(6+2*i , 1)] , [A(4+2*i , 2) A(6+2*i , 2)] ,'Color',Color,'LineWidth',Linewidth)   % H Atome überspringen
        plot( [A(4+2*i , 1) A(5+2*i , 1)] , [A(4+2*i , 2) A(5+2*i , 2)] ,'Color',Color,'LineWidth',Linewidth)   % Verbindungen zu H Atomen
    end

    plot( A(12:14,1) , A(12:14,2),'Color',Color,'LineWidth',Linewidth)                                          % die letzten normalen Verbindungen
    plot( [A(13 , 1) A(1 , 1)] , [A(13 , 2) A(1 , 2)] ,'Color',Color,'LineWidth',Linewidth)                     % Verbindung von Atom 13 zu Atom 1
    plot( [A(3 , 1) A(12 , 1)] , [A(3 , 2) A(12 , 2)] ,'Color',Color,'LineWidth',Linewidth)                     % Verbindung von Atom 3 zu Atom 12
    plot( [A(14 , 1) B(2 , 1)] , [A(14 , 2) B(2 , 2)],'Color',Color,'LineWidth',Linewidth )                     % Verbindung zur nächsten unit cell
    
    A = (R*A')';            % unit cell weiterdrehen
end

plot([0 0] , C(1,:) ,'Color',Color,'LineWidth',Linewidth) 
plot( C(15,:) , [0 0] ,'Color',Color,'LineWidth',Linewidth) 
plot( [0 0]' , C(29,:)' ,'Color',Color,'LineWidth',Linewidth) 
plot( C(43,:) ,  [0 0],'Color',Color,'LineWidth',Linewidth) 


%% Atome dazumalen:

[x,y,z]=sphere(50);

% hold on




 for i=0:3            % Rest plotten
    
  
    
    surf( rN*x-C(1+i*14,1), rN*y-C(1+i*14,2), rN*z+renorm, 'FaceColor', FarbeN, 'Linestyle', 'none' );
    surf( rC*x-C(2+i*14,1), rC*y-C(2+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rC*x-C(3+i*14,1), rC*y-C(3+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rC*x-C(4+i*14,1), rC*y-C(4+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rH*x-C(5+i*14,1), rH*y-C(5+i*14,2), rH*z+renorm, 'FaceColor', FarbeH, 'Linestyle', 'none' );
    surf( rC*x-C(6+i*14,1), rC*y-C(6+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rH*x-C(7+i*14,1), rH*y-C(7+i*14,2), rH*z+renorm, 'FaceColor', FarbeH, 'Linestyle', 'none' );
    surf( rC*x-C(8+i*14,1), rC*y-C(8+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rH*x-C(9+i*14,1), rH*y-C(9+i*14,2), rH*z+renorm, 'FaceColor', FarbeH, 'Linestyle', 'none' );
    surf( rC*x-C(10+i*14,1), rC*y-C(10+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rH*x-C(11+i*14,1), rH*y-C(11+i*14,2), rH*z+renorm, 'FaceColor', FarbeH, 'Linestyle', 'none' );
    surf( rC*x-C(12+i*14,1), rC*y-C(12+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rC*x-C(13+i*14,1), rC*y-C(13+i*14,2), rC*z+renorm, 'FaceColor', FarbeC, 'Linestyle', 'none' );
    surf( rN*x-C(14+i*14,1), rN*y-C(14+i*14,2), rN*z+renorm, 'FaceColor', FarbeN, 'Linestyle', 'none' );
    
   
 end
 
 surf( rCu*x, rCu*y, rCu*z+renorm, 'FaceColor', FarbeCu, 'Linestyle', 'none' );
 
%  lighting gouraud
%  camlight
%  light
%  
%  set(newax,'hittest','off')
%  hlink = linkprop([oldax newax],{'CameraPosition','CameraUpVector','PlotBoxAspectRatio','CameraTarget','DataAspectRatio','PlotBoxAspectRatio','position'});
% key = 'graphics_linkprop';
% setappdata(newax,key,hlink);
 
% lighting gouraud
% light
% light   % oder camlight 
% axis equal
% axis tight
% axis([-dmax dmax -dmax dmax])
% grid off
% xlabel('Angstrøm')
% ylabel('Angstrøm')
% axis off






















































