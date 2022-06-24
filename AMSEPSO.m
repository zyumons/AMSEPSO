function [gbest,gbestval,fitcount]= AMSEPSO(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
rand('state',sum(100*clock));
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
xmin=VRmin;
xmax=VRmax;
Max_FES=10000*D;
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;

fitcount=ps;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);
for i=1:ps
    result(i)=feval(fhd,pos(i,:)',varargin{:});
end

count1=0;count2=0;count3=0;
iwt=0.9-(1:me).*(0.5./me);

vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
lbest=pos;
pbestval=result; %initialize the pbest and the pbest's fitness value
lbestval=result;
[gbestval,gbestid]=min(pbestval);

gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
threshold=0.3*ps;
gbest_stop=0;
cycle=5;
it=1;
estate=1;

while it<=me
    gbest_changed=0;
    
    %singer map
    if it==1
        y=rand;
    end
    mu=0.9+0.18*rand;
    y=mu*(7.86*y-23.31*y.^2+28.75*y.^3-13.30*y.^4);
    if y<0
        y=0;
    end
    c1(it)=2.5-y.*it/me;
    c2(it)=0.5+y.*it/me;
    
    if gbest_stop>=cycle  ||  it==1
        
        [~,gbestid]=min(pbestval);
        gbestr=pbest(gbestid,:);
        n_gbest=0;
        n_gworst=0;
        [~,gworstid]=max(pbestval);
        gworst=pbest(gworstid,:);
        
        for i=1:ps
            if i~=gbestid && i~=gworstid
                
                d1=var(pbest(i,:)-gbestr)/var(gbestr);
                d2=var(pbest(i,:)-gworst)/var(gworst);
                
                if d1>d2
                    n_gworst=n_gworst+1;
                else
                    n_gbest=n_gbest+1;
                end
                
            end
        end
        
        z=n_gworst-n_gbest;
        if z>threshold
            estate=1;
        elseif z<-threshold
            estate=3;
        elseif z>=-threshold && z<=threshold
            estate=2;
        end
        
        if gbest_stop>=cycle
            pp=gbest;
            k=floor(rand*D);
            if k==0
                k=1;
            end
            pp(k)=gbest(k)+(xmax-xmin)*randn;
            if pp(k)>xmax
                pp(k)=xmax;
            elseif pp(k)<xmin
                pp(k)=xmin;
            end
            fitpp=feval(fhd,pp',varargin{:});
            fitcount=fitcount+1;
            if fitpp<gbestval
                gbest=pp;
                gbestval=fitpp;
                gbest_stop=0;
            end
        end
    end
    
    if estate==1
        count1=count1+1;
        p1=pbest;
        p2=lbest;
        
    elseif estate==3
        
        count3=count3+1;
        mbest=sum(pbest)/ps;
        mbestrep=repmat(mbest,ps,1);
        
        [~,index]=sort(pbestval);
        
        for kk=1:ps
            
            a = randi([1,e_ps],1,1);
            b = randi([1,e_ps],1,1);
            c = 0;
            while a == b
                b = randi([1,e_ps],1,1);
            end
            if pbestval(index(a)) < pbestval(index(b))
                cpbest(kk,:) = pbest(index(a),:);
                c = index(a);
            else
                cpbest(kk,:) = pbest(index(b),:);
                c = index(b);
            end
            
            if pbestval(kk) < pbestval(c)
                spbest(kk,:) = pbest(kk,:);
            else
                spbest(kk,:) = cpbest(kk,:);
            end
        end
        
        
        p1=spbest;
        p2=mbestrep;
        
    else
        count2=count2+1;
        
        mbest=sum(pbest)/ps;
        mbestrep=repmat(mbest,ps,1);
        
        p1=pbest;
        p2=mbestrep;
        
    end
    
    for i=1:ps
        aa(i,:)=c1(it).*rand(1,D).*(p1(i,:)-pos(i,:))+c2(it).*rand(1,D).*(p2(i,:)-pos(i,:));
        vel(i,:)=iwt(it).*vel(i,:)+aa(i,:);
        vel(i,:)=(vel(i,:)>mv).*mv+(vel(i,:)<=mv).*vel(i,:);
        vel(i,:)=(vel(i,:)<(-mv)).*(-mv)+(vel(i,:)>=(-mv)).*vel(i,:);
        pos(i,:)=pos(i,:)+vel(i,:);
        pos(i,:)=((pos(i,:)>=VRmin(1,:))&(pos(i,:)<=VRmax(1,:))).*pos(i,:)...
            +(pos(i,:)<VRmin(1,:)).*(VRmin(1,:)+0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D))+(pos(i,:)>VRmax(1,:)).*(VRmax(1,:)-0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D));
        
        if (sum(pos(i,:)>VRmax(i,:))+sum(pos(i,:)<VRmin(i,:)))==0
            result(i)=feval(fhd,pos(i,:)',varargin{:});
            fitcount=fitcount+1;
            if fitcount>=Max_FES
                break;
            end
            tmp=(pbestval(i)<result(i));
            temp=repmat(tmp,1,D);
            pbest(i,:)=temp.*pbest(i,:)+(1-temp).*pos(i,:);
            pbestval(i)=tmp.*pbestval(i)+(1-tmp).*result(i);%update the pbest
            if pbestval(i)<gbestval
                gbest=pbest(i,:);
                gbestval=pbestval(i);
                gbest_changed=1;
            end
        end
    end
    
    for i=1:ps
        if i==1
            if pbestval(ps)<pbestval(i+1)
                lbest(i,:)=pbest(ps,:);
                lbestval(i)=pbestval(ps);
            else
                lbest(i,:)=pbest(i+1,:);
                lbestval(i)=pbestval(i+1);
            end
        elseif i==ps
            if pbestval(i-1)<pbestval(1)
                lbest(i,:)=pbest(i-1,:);
                lbestval(i)=pbestval(i-1);
            else
                lbest(i,:)=pbest(1,:);
                lbestval(i)=pbestval(1);
            end
        else
            if pbestval(i-1)<pbestval(i+1)
                lbest(i,:)=pbest(i-1,:);
                lbestval(i)=pbestval(i-1);
            else
                lbest(i,:)=pbest(i+1,:);
                lbestval(i)=pbestval(i+1);
            end
        end
    end
    
    if gbest_changed==1  % 更新 gbest 停滞代数
        gbest_stop=0;
    else
        gbest_stop=gbest_stop+1;
    end
    if fitcount>=Max_FES
        break;
    end
    if (it==me)&&(fitcount<Max_FES)
        it=it-1;
    end
    it=it+1;
end