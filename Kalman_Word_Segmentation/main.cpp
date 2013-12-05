//
//  main.cpp
//  Kalman_on_Segment
//
//  Created by Zhu Weimeng on 12-10-30.
//  Copyright (c) 2012年 Zhu Weimeng. All rights reserved.
//

#include <algorithm>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <map>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <string>
#define ch unsigned char
#define sqr(x) ((x)*(x))
//#define DEBUG
using namespace std;
//#pragma comment(linker, "/STACK:500000000 ")
int n,sum1,sum2,f[256][256];
ch * s;
ch * s0;
ch * s2;
#define w 10

class cmp {
public:
	bool operator()(const ch* c1,const ch* c2) const {
		for (int i=0;i<4;i++)
			if (c1[i]!=c2[i]) return bool(c1[i]<c2[i]);
		return 0;
	}
};


double * mi,*dts,*md;
map<ch*,int,cmp> f2;
map<ch*,double,cmp>::iterator iter;

long len_arti;


void preinit() {
	FILE * article;
	article=fopen("data.txt","r");
	fseek(article,0L,SEEK_END);
	len_arti = ftell(article);
	len_arti*=1.5;
	//len=(long)floor(len/1.5);
	s= new ch[len_arti];
	s0= new ch[len_arti];
	mi= new double[len_arti];
	dts = new double[len_arti];
	md = new double[len_arti];
	s2= new ch[len_arti];
	fclose(article);
	printf("len_arti: %ld\n",len_arti);
}



void prein() {  //对于输入进行预处理，进行虚词处理
	map<string,int> my_map;
	map<string,int>::iterator my_pointer;
	unsigned char a[3];
	char p[3];
	bool flag=false;
	int n=0;
    FILE * dict_txt;
    dict_txt=fopen("dict.txt","r");
	//freopen("dict.txt","r",stdin);
	while (fscanf(dict_txt,"%s",p) != EOF)
		my_map.insert(pair<string,int>(p,0));
	//freopen("旧唐书.txt","r",stdin);
	a[2] = '\0';
	p[2] = '\0';
    FILE * article;
    article=fopen("data.txt","r");
	while (fscanf(article,"%c",&a[0]) != EOF) {
		while (s0[n]) n++;
		if (a[0] == '\n') {
			memcpy(s0+n,"\n",1);
			flag =false;
			continue;
		}
		if (a[0] < 127) {
			if (a[0] != ' ') {
				memcpy(s0+n,a,1);
				flag=true;
			}
			continue;
		}
		fscanf(article,"%c",&a[1]);
		p[0] = a[0];
		p[1] = a[1];
		if (a[0] >= 0xb0 && a[1] >= 0xa1)
		{
			my_pointer = my_map.find(p);
			if (my_pointer == my_map.end()) {
				memcpy(s0+n,p,sizeof(p));
				flag=false;
			}
			else {
				if (!flag) {
					flag = true;
					memcpy(s0+n,"/",1);
				}
				while (s0[n]) n++;
				memcpy(s0+n,p,sizeof(p));
				while (s0[n]) n++;
				memcpy(s0+n,"/",1);
				flag = true;
			}
		} else {
			memcpy(s0+n,p,sizeof(p));
			flag=true;
		}
	}
    fclose(dict_txt);
    fclose(article);
}

bool check(int i) {  //检查s中第i个字节开始的字符是否是汉字
	if (i+1>=n) return 0;
	return (s[i]>=0xB0 && s[i+1]>=0xA1 && s[i]<=0xF7 && s[i+1]<=0xFE)
    ||(s[i]>=0x81 && s[i]<=0xA0 && ((s[i+1]>=0x40 && s[i+1]<=0x7E) || (s[i+1]>=0x80 && s[i+1]<=0xFE)))
    ||(s[i]>=0xAA && s[i]<=0xFE && ((s[i+1]>=0x40 && s[i+1]<=0x7E) || (s[i+1]>=0x80 && s[i+1]<=0xA0)));
}

bool check(ch* s,int i,int n){
	if (i+1>=n) return 0;
	return (s[i]>=0xB0 && s[i+1]>=0xA1 && s[i]<=0xF7 && s[i+1]<=0xFE)
    ||(s[i]>=0x81 && s[i]<=0xA0 && ((s[i+1]>=0x40 && s[i+1]<=0x7E) || (s[i+1]>=0x80 && s[i+1]<=0xFE)))
    ||(s[i]>=0xAA && s[i]<=0xFE && ((s[i+1]>=0x40 && s[i+1]<=0x7E) || (s[i+1]>=0x80 && s[i+1]<=0xA0)));
}

void init() {	//读入整个旧唐书文本，存入s中，为了计算方便将单字节字符添一个0补充成双字节
	int p,j;
	int i=0;
	s[0]=10;s[1]=0;n=2;
	while (p = s0[i++]) {
		if (p<128) {
			if (p==47) {
				s[n]=s0[i];
				s[n+1]=s0[i+1];
				j=n;n+=10000;
				if (!check(j)) {
					n-=10000;
					continue;
				}
				n-=10000;
			}
			s[n++]=p;
			s[n++]=0;
		}else {
			s[n++]=p;
			p=s0[i++];
			s[n++]=p;
		}
	}
	s[n++]=0;s[n++]=0;
}

void calcf1() { //计算所有单个汉字出现的频数（为了计算t-test字符也进行了统计）
	for (int i=0;i<n;i+=2) {
		f[s[i]][s[i+1]]++;
		sum1++;
	}
}

void adjust(int k) {  //调整在顿号所夹的词的互信息值
	char s[3];
	int i,flag;
	unsigned char cha[3],a[100];
	//freopen("旧唐书.txt","r",stdin);
	//freopen("dict_d.txt","w",stdout);
    FILE * article;
    article=fopen("data.txt","r");
    
    FILE * dict_d_txt;
    dict_d_txt=fopen("dict_d.txt", "w");
    
	cha[2] = '\0';s[2] = '\0';
	while (fscanf(article,"%c",&cha[0]) != EOF) {
		if (cha[0] < 127) continue;
		fscanf(article,"%c",&cha[1]);
		if (cha[0] >= 0xb0 && cha[1] >= 0xa1) continue;
		s[0] = cha[0];
		s[1] = cha[1];
		if (strcmp(s,"、") == 0) {
			flag = 0;
			for (i = 0;i <= 4;i++) {
				fscanf(article,"%c",&cha[0]);
				if (cha[0] < 127) break;
				fscanf(article,"%c",&cha[1]);
				a[2 * i] = cha[0];
				a[2 * i + 1] = cha[1];
				s[0] = cha[0];
				s[1] = cha[1];
				if (strcmp(s,"、") == 0) {
					flag = 1;
					break;
				};
				if (!(cha[0] >= 0xb0   &&  cha[1] >= 0xa1) )
					break;
			}
			a[i * 2] = '\0';
			if ((flag == 1) && (i > 1)) {
				for (int j=0;a[j];j+=2) {
					ch *p=new ch(5);
					memcpy(p,a+j,4);
					p[4]=0;
					if (f2[p]*k>0) {
						if (k>0) f2[p]*=-32;
						else f2[p]/=-32;
					}
				}
				if (k>0) fprintf(dict_d_txt,"%s\n",a);
			}
		}
	}
    fclose(article);
    fclose(dict_d_txt);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//所有词典，公用一个map，二元组，键是词，值是每个词典的lambda，可以自定义
map<string,double> dicts;
map<string,double>::iterator dicts_pointer;
//几个规则：
//0、lambda不要大于0.5，要不然规则3的加值就麻烦了??
//1、如果是单字的词应该倾向于切开，把lambda变成负的再添加进入词典，不过也可以让-lambda不要太小??
//2、如果是双字的词直接加入词典，词典中的值即是lambda
//3、如果是多于两个字的词，则给前两个字建一个词条，然后整数部分是这个词的长度，小数部分是0
//		在完整词的词条项目上值是lambda（或者随着词的位数增大越来越大，比如lambda+1/2lambda+1/4lambda+...）
//4、如果一个词既是一个词又是另一个词的前缀，则整数部分是前缀所在词的长度，小数部分是lambda
//5、如果一个词是n多词的前缀，则整数部分是最长的前缀词所在的长度
//6、前缀是单字词的与多字词相同只不过是负数。
//7、两个词典里出现了同一个词，取lambda大的加入词典
//8、不考虑词典里的词不是汉字的情况??
//把词典读到s2里边
void adddict(char* dict,double lambda)	//添加词典
{
	memset(s2,0,len_arti);
	//freopen(dict,"r",stdin);
    FILE * tmpdict = fopen(dict, "r");
	ch a[3];
	a[2] = '\0';
	char b[3];
	b[2] = '\0';
	int n = 0;
	while (fscanf(tmpdict,"%c",&a[0]) != EOF)		//一个词一行
	{
		while (s2[n]) n++;
		if (a[0] == '\n')
		{
			memcpy(s2+n,"\n",1);
			//flag =false;
			continue;
		}
		fscanf(tmpdict,"%c",&a[1]);
		if (check(a,0,2))
		{
			b[0] = a[0];
			b[1] = a[1];
			memcpy(s2+n,b,sizeof(b));
		}
	}
	while (s2[n]) n++;
	int len = 0;
	string p = "";
	for (int i=0;i<=n;i+=2)
	{
		if (check(s2,i,n))						//读取到一整个词之后再处理
		{
			len++;
			p += string(1,s2[i])+string(1,s2[i+1]);
		}
		else
		{
			if(len == 0)
				continue;
			else if(len == 1)				//1、如果是单字的词应该倾向于切开，把lambda变成负的再添加进入词典
			{
				dicts_pointer = dicts.find(p);
				if(dicts_pointer == dicts.end())
					dicts.insert(pair<string,double>(p,0-lambda));
				else
				{
					if(dicts_pointer->second > 0-lambda && dicts_pointer->second > -1)	//之前加进去的是单个词
						dicts_pointer->second = 0-lambda;
					else if(dicts_pointer->second <= -1)														//此单字词是其他词的前缀，比较一下小数部分的大小
					{
						if(ceil(dicts_pointer->second) - dicts_pointer->second < lambda)
							dicts_pointer->second = ceil(dicts_pointer->second)-lambda;
					}
				}
			}
			else if(len == 2)				//2、如果是双字的词
			{
				dicts_pointer = dicts.find(p);
				if(dicts_pointer == dicts.end())								//以前没出现过直接加入词典，词典中的值即是lambda
					dicts.insert(pair<string,double>(p,lambda));
				else
				{
					if(dicts_pointer->second >= 1)							//此双字词是其他词的前缀，比较小数部分大小
					{
						if(dicts_pointer->second - floor(dicts_pointer->second) < lambda)
							dicts_pointer->second = floor(dicts_pointer->second)+lambda;
					}
					if(dicts_pointer->second < lambda)
						dicts_pointer->second = lambda;
				}
			}
			else
			{
				string p1(p,0,2);		//取单字前缀
				dicts_pointer = dicts.find(p1);
				if(dicts_pointer != dicts.end())			//如果确实是个单字词（如果这个词不是单字词就不用管它了）
				{
					if(dicts_pointer->second > -1)		//如果他不是别的词的前缀
						dicts_pointer->second -= len;
					else												//这个单字词已经是其他词的前缀了，比较这两个词哪个长，把整数部分改成最长的词
					{
						if(ceil(dicts_pointer->second) > 0-len)
							dicts_pointer->second = dicts_pointer->second - ceil(dicts_pointer->second) - len;
					}
				}
				string p2(p,0,4);		//取两字前缀
				dicts_pointer = dicts.find(p2);
				if(dicts_pointer != dicts.end())			//如果确实是个两字词，或者是别的词的前缀
				{
					if(dicts_pointer->second < 1)		//如果他不是别的词的前缀
						dicts_pointer->second += len;
					else												//他已经是别的词的前缀了，比较一下两个词哪个长，把整数部分改成最长的词
					{
						if(floor(dicts_pointer->second) < len)
							dicts_pointer->second = dicts_pointer->second - floor(dicts_pointer->second)+len;
					}
				}
				else													//这个词不是一个两字词，把前两个字作为前缀插入
					dicts.insert(pair<string,double>(p2,len));
                
				//前缀处理好了，把这个词整个加进去
				double longlambda = lambda;
				for(int i = 3;i <= len;i++)
					longlambda += (1/pow(2.0,i-2))*lambda;
				dicts_pointer = dicts.find(p);
				if(dicts_pointer != dicts.end())			//如果以前插入过这个词
				{
					if(dicts_pointer->second < longlambda)//取大的lambda
						dicts_pointer->second = longlambda;
				}
				else
					dicts.insert(pair<string,double>(p,longlambda));
			}
			len = 0;
			p = "";
			i--;
		}
	}
    fclose(tmpdict);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void calcmi() { //计算每相邻两个汉字的互信息，存为mi数组，暂时没有加入可信度的调整
	for (int i=0;i<n;i+=2)
		if (check(i) && check(i+2)) //计算二元组的总数
			sum2++;
	for (int i=0;i<n;i+=2)
		if (check(i) && check(i+2)) {
			ch *c=(ch*)malloc(5);
			for (int j=0;j<4;j++) c[j]=s[i+j]; //计算每个二元组的个数
			c[4]=0;
			f2[c]++;
		}
    adjust(1);
    for (int i=0;i<n;i+=2)
        if (check(i) && check(i+2)) { //计算互信息mi值
            ch *c=(ch*)malloc(5);
            for (int j=0;j<4;j++) c[j]=s[i+j];
            c[4]=0;
            double d=(double)abs(f2[c])*sum1/f[s[i]][s[i+1]]/f[s[i+2]][s[i+3]]/sum2*sum1;
            d=log((double)d)/log((double)2);
            mi[i]=d;
        }
    adjust(-1);
}

void calcdts() { //计算相邻两个汉字的dts值，为左右两个三元组的t-test差
	for (int i=0;i<n-4;i+=2)
		if (check(i+2) && check(i+4)) {
			ch *c=(ch*)malloc(15);
			for (int j=0;j<4;j++) c[j]=s[i+j];
			for (int j=0;j<4;j++) c[5+j]=s[i+j+2];
			for (int j=0;j<4;j++) c[10+j]=s[i+j+4];
			c[4]=0;c[9]=0;c[14]=0;
			double p[3];
			for (int j=0;j<3;j++) p[j]=f2[c+j*5];
			int f1[4];
			for (int j=0;j<4;j++) f1[j]=f[s[i+j*2]][s[i+j*2+1]];
			double t1,t2;
			t1=(p[0]/f1[0]-p[1]/f1[1])/sqrt(p[0]/f1[0]/f1[0]+p[1]/f1[1]/f1[1]);
			t2=(p[1]/f1[1]-p[2]/f1[2])/sqrt(p[1]/f1[1]/f1[1]+p[2]/f1[2]/f1[2]);
			dts[i+2]=t2-t1;
		}
}

void calcmd(double t1,double t2,double lambda){  //对mi值和dts值进行线性叠加，t1和t2分别为mi和dts的阀值
	double avr_mi=0,avr_dts=0,sigma_mi=0,sigma_dts=0;
	for(int i=2;i<n-4;i+=2)  //计算mi和dts的平均值
		if(check(i)&&check(i+2)){
			avr_mi+=mi[i]/(n/2-3);
			avr_dts+=dts[i]/(n/2-3);
		}
    for(int i=2;i<n-4;i+=2)  //计算mi和dts的方差
        if(check(i)&&check(i+2)){
            sigma_mi+=sqr(mi[i]-avr_mi)/(n/2-3);
            sigma_dts+=sqr(dts[i]-avr_dts)/(n/2-3);
        }
    sigma_mi=sqrt(sigma_mi);  //计算标准差
    sigma_dts=sqrt(sigma_dts);
    for(int i=2;i<n-4;i+=2)
        if(check(i)&&check(i+2))
            md[i]=(mi[i]-t1)/sigma_mi+lambda*(dts[i]-t2)/sigma_dts;  //进行线性叠加
}

bool vcmp(pair<double,ch*> a,pair<double,ch*> b) { //调试时按照互信息与dts值输出时的比较函数
	return a.first>b.first;
}

void showmi() { //按照互信息值排序输出二元字组(调试用)
	FILE *MI=fopen("swudmi.out","w");
	vector<pair<double,ch*> > t;
	for (int i=2;i<n-4;i++){
		ch *c=new ch[5];
		c[0]=s[i];c[1]=s[i+1];c[2]=s[i+2];c[3]=s[i+3];c[4]=0;
		t.push_back(make_pair(mi[i],c));
	}
	sort(t.begin(),t.end(),vcmp);
	for (int i=0;i<t.size();i++)
		fprintf(MI,"%s,%.6lf\n",t[i].second,t[i].first);
    fclose(MI);
}

void showdts() { //按照dts值排序输出二元字组（调试用）
	FILE *DTS=fopen("swuddts.out","w");
	vector<pair<double,ch*> > t;
	for (int i=2;i<n-4;i++){
		ch *c=new ch[5];
		c[0]=s[i];c[1]=s[i+1];c[2]=s[i+2];c[3]=s[i+3];c[4]=0;
		t.push_back(make_pair(dts[i],c));
	}
	sort(t.begin(),t.end(),vcmp);
	for (int i=0;i<t.size();i++)
		fprintf(DTS,"%s,%.6lf\n",t[i].second,t[i].first);
    fclose(DTS);
}

void showmd() {  //按照md值排序输出二元字组（调试用）
	FILE *MD=fopen("swudmd.out","w");
	vector<pair<double,ch*> > t;
	for (int i=2;i<n-4;i++){
		ch *c=new ch[5];
		c[0]=s[i];c[1]=s[i+1];c[2]=s[i+2];c[3]=s[i+3];c[4]=0;
		t.push_back(make_pair(md[i],c));
	}
	sort(t.begin(),t.end(),vcmp);
	for (int i=0;i<t.size();i++)
		fprintf(MD,"%s,%.6lf\n",t[i].second,t[i].first);
    fclose(MD);
}

bool checknum(ch *a,ch number[][3]){  //检查一个字是否为中文数字
	for(int i=0;i<14;i++)
		if(a[0]==number[i][0]&&a[1]==number[i][1]) return true;
	return false;
}

bool compare(ch *a,char tmp[2]){  //将两个中文汉字进行比较
	ch word[2];
	memcpy(word,tmp,sizeof word);
	if(a[0]==word[0]&&a[1]==word[1]) return true;
	return false;
}

void printDict()
{
	dicts_pointer = dicts.begin();
	//freopen("dictionary.txt","w",stdout);
    FILE * prdict = fopen("dictionary.txt", "w");
	int num = 1;
	while(dicts_pointer != dicts.end())
	{
		fprintf(prdict,"%d %s %f\n",num,dicts_pointer->first.c_str(),dicts_pointer->second);
		dicts_pointer++;
		num++;
	}
    fclose(prdict);
}

/*********************By Zhu Weimeng START**********************************************************/

double minest(double a, double b, double c) {
    if (a<b && a<c) {
        return a;
    }
    else if(b<c) {
        return b;
    }
    else
        return c;
}

double maxest(double a, double b, double c) {
    if (a>b && a>c) {
        return a;
    }
    else if(b>c) {
        return b;
    }
    else
        return c;
}

void begin_kalman(ch * s, double * val, double t) {
    FILE* testfile = fopen("error_output.txt", "a");
    FILE* testinput = fopen("test_input.txt", "r");
    FILE* modeloutput = fopen("model_output.txt","a");
    
    //only allocate 1/100 of the max text character count to save space
    double *target_md;
    int *target_address;
    double *upper_context;
    double *latter_context;
    double max_peak=-100.0, min_bottom=100.0;
    
    //use to count error stats
    int error_kalman=0;
    int error_dts=0;
    int error_md=0;
    int error_mi=0;
    
    //get target word
    ch* target = new ch[10];
    while(fscanf(testinput,"%s",target)!=EOF) {
        //while(scanf("%s",target) && target[0]!='e') {
        
        target_md=new double[50000];
        target_address=new int[50000];
        upper_context=new double[50000];
        latter_context=new double[50000];
        max_peak=0.0-100.0;
        min_bottom=100.0;
        
        //use to count error stats
        error_kalman=0;
        error_dts=0;
        error_md=0;
        error_mi=0;
        
        printf("target word: %s\n",target);
        //scanf("%s",target);
        fprintf(testfile,"\ntarget word: %s\n",target);
        fprintf(modeloutput,"\ntarget word: %s\n",target);
        fprintf(modeloutput," X_pre K_Gain X_post   Z\n");
        
        //printf("init peak:%1.4lf bottom:%1.4lf\n",max_peak,min_bottom);
        
        //find target_md distribution and context md value&max&min
        int target_num=0;
        for (int i=2; i<n-2; i+=2) {
            if (target[0]==s[i-2] && target[1]==s[i-1] && target[2]==s[i] && target[3]==s[i+1]) {
                target_md[target_num]=val[i-2];
                target_address[target_num]=i;
                if(s[i-4]<127)
                    upper_context[target_num]=val[i-6];
                else {
                    upper_context[target_num]=val[i-4];
                }
                if (upper_context[target_num]>max_peak && check(i-4) && upper_context[target_num]<100.0) {
                    max_peak=upper_context[target_num];
                }
                if (latter_context[target_num]>max_peak && check(i+2) && latter_context[target_num]<100.0) {
                    max_peak=latter_context[target_num];
                }
                if(s[i+2]<127)
                    latter_context[target_num]=val[i+2];
                else {
                    latter_context[target_num]=val[i];
                }
                if (latter_context[target_num]<min_bottom && check(i+2) && latter_context[target_num]>-100.0) {
                    min_bottom=latter_context[target_num];
                }
                if (upper_context[target_num]<min_bottom && check(i-4)  && upper_context[target_num]>-100.0) {
                    min_bottom=upper_context[target_num];
                }
                
                //printf("uc:%lf lc:%lf\n",upper_context[target_num],latter_context[target_num]);
                
                target_num++;
            }
        }
        
        //printf("pre peak:%1.4lf bottom:%1.4lf\n",max_peak,min_bottom);
        
        if(max_peak<(t+0.5))
            max_peak=(t+0.5);
        if(min_bottom>(t-0.5))
            min_bottom=(t-0.5);
        
        //printf("post peak:%1.4lf bottom:%1.4lf\n",max_peak,min_bottom);
        
        printf("target_num: %d\n",target_num);
        
        //calc Q of the hypothesis W~N(0,Q)
        double md_bar=0.0;
        for (int i=0; i<target_num; i++) {
            md_bar+=target_md[i];
        }
        md_bar/=target_num;
        double Q = 0.0;
        for (int i=0; i<target_num; i++) {
            Q+=(target_md[i]-md_bar)*(target_md[i]-md_bar);
        }
        Q/=target_num;
        
        //start kalman filter
        double X=md_bar;
        double P=Q;
        double Z;
        double R;
        double K;
        double W=0.0;
        
        fprintf(testfile, "kalman   dts    mi   ori_md\n");
        
        for (int i=0; i<target_num; i++) {
            //do predict
            /**************CORE**************/
            X=X+W;
            P=P+Q;
            /**************CORE*END**********/
            
            fprintf(modeloutput, "%1.4lf ",X);
            
            printf("%d: ",i+1);
            
            //show original seg result
            bool seg_mark=false;
            for (int j=target_address[i]-40; j<target_address[i]+40; j+=2) {
                if (j<0) {
                    continue;
                }
                if (j>=n) {
                    break;
                }
                if(j==target_address[i]-2)
                    printf(" ");
                if(j==target_address[i]+2)
                    printf(" ");
                
                if (s[j]<127) {
                    printf("%c",s[j]);
                }
                else if (!check(j)) {
                    printf("%c%c",s[j],s[j+1]);
                }
                else {
                    if (j==target_address[i]) {
                        if(check(j-4)&&check(j+2)&&X>=upper_context[i]&&X>=latter_context[i]&&X<t-0.5) {
                            printf("/");  //波峰情况的判定
                            seg_mark=true;
                        }
                        else if(check(j-4)&&check(j+2)&&X<=upper_context[i]&&X<=latter_context[i]&&X<=t+0.5) {
                            printf("/");  //波谷情况的判定
                            seg_mark=true;
                        }
                        else if(X<=t) {
                            printf("/");  //一般情况判定
                            seg_mark=true;
                        }
                        else {
                            printf("-");
                            seg_mark=false;
                        }
                    }
                    printf("%c%c",s[j],s[j+1]);
                }
            }
            
            //original dts result
            bool dts_seg=false;
            if (dts[target_address[i]-4]>=dts[target_address[i]-2] && dts[target_address[i]]>=dts[target_address[i]-2]) {
                dts_seg=true;
            }
            
            //original mi result
            bool mi_seg=false;
            if (mi[target_address[i]-4]>=mi[target_address[i]-2] && mi[target_address[i]]>=mi[target_address[i]-2]) {
                mi_seg=true;
            }
            
            //original md result
            bool md_seg=false;
            if(check(target_address[i]-4)&&check(target_address[i]+2)&&val[target_address[i]-2]>=val[target_address[i]-4]&&val[target_address[i]-2]>=val[target_address[i]]&&val[target_address[i]-2]<t-0.5) {
                md_seg=true;  //波峰情况的判定
            }
            else if(check(target_address[i]-4)&&check(target_address[i]+2)&&val[target_address[i]-2]<=val[target_address[i]-4]&&val[target_address[i]-2]<=val[target_address[i]]&&val[target_address[i]-2]<=t+0.5) {
                md_seg=true;  //波谷情况的判定
            }
            else if(val[target_address[i]-2]<=t) {
                md_seg=true;  //一般情况判定
            }
            
            printf("\n\n:");
            //get observation
            char oper[2];
            scanf("%s",oper);
            if (oper[0]=='c') {
                //correct seg
                Z=X;
                //split correct
                if (seg_mark) {
                    if (X>=upper_context[i]&&X>=latter_context[i]&&X<t-0.5) {
                        R=(X-(t-0.5))*(X-(t-0.5))*4/12;
                    }
                    else if (X<=upper_context[i]&&X<=latter_context[i]&&X<=t+0.5) {
                        R=(X-minest(upper_context[i],latter_context[i],t+0.5))*(X-minest(upper_context[i],latter_context[i],t+0.5))*4/12;
                    }
                    else if (X<=t) {
                        R=(X-t)*(X-t)*4/12;
                    }
                }
                //combine correct
                else {
                    R=(X-t)*(X-t)*4/12;
                }
                
                //other method count
                if (seg_mark!=dts_seg) {
                    error_dts++;
                }
                if (seg_mark!=mi_seg) {
                    error_mi++;
                }
                if (seg_mark!=md_seg) {
                    error_md++;
                }
            }
            if (oper[0]=='n') {
                //incorrect seg
                if (seg_mark) {
                    double Z_min=t;
                    if (X<=upper_context[i]&&X<=latter_context[i]&&X<=t+0.5) {
                        Z_min=t+0.5;
                    }
                    Z=(Z_min+max_peak)/2.0;
                    R=(max_peak-Z_min)*(max_peak-Z_min)/12.0;
                }
                else {
                    double Z_max;
                    if ((t+0.5)<=upper_context[i]&&(t+0.5)<=latter_context[i]) {
                        Z_max=t+0.5;
                    }
                    else if((t-0.5)>=upper_context[i]&&(t-0.5)>=latter_context[i]) {
                        Z_max=t-0.5;
                    }
                    else {
                        Z_max=t;
                    }
                    Z=(min_bottom+Z_max)/2.0;
                    R=(Z_max-min_bottom)*(Z_max-min_bottom)/12.0;
                }
                
                //kalman error count
                error_kalman++;
                //other method error count
                if (seg_mark==dts_seg) {
                    error_dts++;
                }
                if (seg_mark==md_seg) {
                    error_md++;
                }
            }
            if(oper[0]=='e') {
                fclose(testinput);
                fclose(testfile);
                fclose(modeloutput);
                return;
            }
            
            //output error result
            fprintf(testfile,"%1.4lf %1.4lf %1.4lf %1.4lf\n",error_kalman/(i+1.0),error_dts/(i+1.0),error_mi/(i+1.0),error_md/(i+1.0));
            
            //reverse predict
            /**************CORE**************/
            K=P/(P+R);
            X=X+K*(Z-X);
            P=(1-K)*P;
            /**************CORE*END**********/
            fprintf(modeloutput, "%1.4lf %1.4lf %1.4lf\n",K,X,Z);
        }
        
        delete target_md;
        delete target_address;
        delete upper_context;
        delete latter_context;
        
        printf("end of %s\n\n", target);
        
    }
    
    fclose(testinput);
    fclose(testfile);
    fclose(modeloutput);
}

/*********************By Zhu Weimeng END************************************************************/

int main() {
    printf("preinit()\n");
	preinit();
	printf("prein()\n");
	prein();
	printf("init()\n");
	init();
	printf("adddict()\n");
	adddict("dict_place.txt",0.5);
	adddict("dict_people.txt",0.5);
	adddict("dict_office.txt",0.5);
	printf("calcf1()\n");
	calcf1();
	printf("calcmi()\n");
	calcmi();
	printf("calcdts()\n");
	calcdts();
	printf("calcmd()\n");
	calcmd(5,0,0.7);
#ifdef DEBUG
	showmi();
	showdts();
	showmd();
#endif
    printf("begin_kalman()\n");
    begin_kalman(s,md,0.3);
	//printDict();
	//print(s,md,0.3);
	return 0;
}

