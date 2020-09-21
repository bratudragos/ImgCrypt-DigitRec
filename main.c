#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
struct imgbmp
{
    unsigned char *infoheader;
    unsigned char *pixels;
    unsigned char *padding;
    unsigned int width, height,size,paddingsize;
};

struct coord{
    long *x,*y,lungime;
    float *corr;
    int *cifra;
};

struct coord2{
    long x,y,lungime;
    float corr;
    int cifra;
};

uint32_t octeti(uint32_t x,int a){
    if(a==0){
        return (x & 255);
    }
    if(a==1){
        x= (x & 65280);
        x=x>>8;
        return x;
    }
    if(a==2)
    {
        x= (x & 16711680);
        x=x>>16;
       return x;
    }
    return 0;
}

uint32_t xorshift32(uint32_t *state)
{
    uint32_t x = state[0];
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    state[0] = x;
    return x;
}

void elibmemorie(struct imgbmp *imglin)
{
    if(imglin->infoheader!=NULL)
        free(imglin ->infoheader);
    if(imglin->pixels!=NULL)
        free(imglin->pixels);
    if(imglin->padding!=NULL)
        free(imglin->padding);
    if(imglin!=NULL)
        free(imglin);
}

void elibmemorie2(struct coord *ferestre)
{
    if(ferestre->x!=NULL)
      free(ferestre->x);
    if(ferestre->y!=NULL)
        free(ferestre->y);
    if(ferestre->corr!=NULL)
        free(ferestre->corr);
    if(ferestre->cifra!=NULL)
        free(ferestre->cifra);
    if(ferestre!=NULL)
        free(ferestre);
}

struct imgbmp* incarca(char *str)
{

    unsigned int i,j,l;
    unsigned long k=0,h=0;
    FILE *f = fopen(str, "rb");
    if (f==NULL)
    {

        printf("eroare la deschidere imagine");
        return 0;
    }
        struct imgbmp *imglin;
        imglin =(struct imgbmp*)malloc(sizeof(struct imgbmp));
        imglin->infoheader=(unsigned char*)malloc(54);
        fread(imglin->infoheader, sizeof(unsigned char), 54, f);
        imglin->width = *(unsigned int *)&(imglin->infoheader)[18];
        imglin->height = *(unsigned int *)&(imglin->infoheader)[22];
        imglin->size=imglin->height*imglin->width*3;
        //imglin->size+=((4-((imglin->width*3)%4))*imglin->height);
        imglin->pixels =(unsigned char*)malloc(imglin->size);
        imglin->paddingsize=4-((imglin->width*3)%4);
        if(imglin->paddingsize==4)
            imglin->paddingsize=0;
        imglin->padding=(unsigned char*)malloc(imglin->paddingsize*imglin->height);
        unsigned char c;
        //printf("\n%d\n",imglin->height);
        for(i=0;i<imglin->height;i++)
        {
            for(j=0;j<imglin->width;j++){
                for(l=0;l<3;l++)
                {fread(&c,1,1,f);
                imglin->pixels[k++]=c;
                }
            }
            for(j=0;j<imglin->paddingsize;j++)
            {
                fread(&c,1,1,f);
                imglin->padding[h++]=c;
            }

        }
        //fread(imglin->pixels, sizeof(unsigned char),imglin->size, f);
        fclose(f);
        return imglin;

}

void salvare(char *str,struct imgbmp *imglin){
  FILE *f=fopen(str,"wb");
  long i,j,l,k=0,aux=0;
  unsigned char c;
  fwrite(imglin->infoheader,sizeof(unsigned char),54,f);
  /*for(i=0;i<imglin->height;i++){
      fwrite(imglin->pixels+(i*imglin->width*3),1,imglin->width*3,f);
      fwrite(imglin->padding+(i*imglin->paddingsize),1,imglin->paddingsize,f);

}*/
    for(i=0;i<imglin->height;i++)
    {
        for(j=0;j<imglin->width;j++){
            for(l=0;l<3;l++)
            {c=imglin->pixels[k];
                fwrite(&c,1,1,f);
                k++;
            }
        }
        for(j=0;j<imglin->paddingsize;j++)
        {
            c=imglin->padding[aux];
            fwrite(&c,1,1,f);
            aux++;
        }

    }
    fclose(f);
}

void criptare(char *init,char *cript,char *cheie){
    FILE *f=fopen(init,"rb");
    FILE *g=fopen(cript,"wb");
    FILE *key=fopen(cheie,"r");
// incarca imaginea
    struct imgbmp *imglin=incarca(init);
//nr aleatoare cu xorshift
    uint32_t *state=(uint32_t*)malloc(sizeof(uint32_t));
    uint32_t aux,sv=0;
    fscanf(key,"%d",&aux);
    //printf("%d ",aux);
    fscanf(key,"%d",&sv);
    //printf("%d",sv);

    //fscanf(key,"%d",&sv);
    //printf("%d",sv);
    state[0]=aux;
    uint32_t *secv=(uint32_t*)malloc(((2*imglin->width*imglin->height))*sizeof(uint32_t));
    int i,j;
    for (i = 0; i < (2*imglin->width*imglin->height)-1; i++)
    {
        uint32_t x=xorshift32(state);
        secv[i]=x;
    }
//determinarea permutarii
    uint32_t *permut=(uint32_t*)malloc(imglin->width*imglin->height*sizeof(uint32_t));
    uint32_t p;
    for(p=0;p<imglin->width*imglin->height;p++)
        permut[p]=p;

    for(i=(imglin->height*imglin->width)-1;i>=0;i--){
      j=secv[i]%(i+1);
      aux=permut[i];
      permut[i]=permut[j];
      permut[j]=aux;
    }

//permutarea pixelilor
    struct imgbmp *imgtemp=(struct imgbmp*)malloc(sizeof(struct imgbmp));
    imgtemp->infoheader=(unsigned char*)malloc(54);
    //strcpy(imgtemp->infoheader,imglin->infoheader);
    for(i=0;i<54;i++){
            imgtemp->infoheader[i]=imglin->infoheader[i];
       }
    imgtemp->width = *(unsigned int *)&(imglin->infoheader)[18];
    imgtemp->height = *(unsigned int *)&(imglin->infoheader)[22];
    imgtemp->size=imgtemp->width*imgtemp->height*3;
    imgtemp->pixels = malloc(imgtemp->size);
    imgtemp->paddingsize=4-((imgtemp->width*3)%4);
    if(imgtemp->paddingsize==4)
        imgtemp->paddingsize=0;
    imgtemp->padding=(unsigned char*)malloc(imgtemp->paddingsize*imgtemp->height);

    for(i=0;i<imgtemp->width*imgtemp->height;i++){
        {imgtemp->pixels[permut[i]*3]=imglin->pixels[i*3];
         imgtemp->pixels[permut[i]*3+1]=imglin->pixels[i*3+1];
         imgtemp->pixels[permut[i]*3+2]=imglin->pixels[i*3+2];
        }
    }

    /*
     * p= pr,pg,pb
     * uint32 x
     * p^x= pr^x2, pg^x1, pb^x0
     * x2=x & 16 711 680
     * x1=x & 65 280
     * x0=x & 255
     *
     */
    // criptarea in sine
   struct imgbmp *imgtemp2=(struct imgbmp*)malloc(sizeof(struct imgbmp));
    imgtemp2->infoheader=(unsigned char*)malloc(54);
    //strcpy(imgtemp2->infoheader,imglin->infoheader);
    for(i=0;i<54;i++){
        imgtemp2->infoheader[i]=imglin->infoheader[i];
    }
    imgtemp2->width = *(unsigned int *)&(imglin->infoheader)[18];
    imgtemp2->height = *(unsigned int *)&(imglin->infoheader)[22];
    imgtemp2->size=imglin->height*imglin->width*3;
    imgtemp2->pixels = malloc(imglin->size);
    imgtemp2->paddingsize=4-((imgtemp2->width*3)%4);
    if(imgtemp2->paddingsize==4)
        imgtemp2->paddingsize=0;
    imgtemp2->padding=(unsigned char*)malloc(imgtemp2->paddingsize*imgtemp2->height);


    //primul pixel
     /*
      * sv0=octeti(sv,0)
      * sv1=octeti(sv,1)
      * sv2=octeti(sv,2)
      * r0=octeti(r,0)
      * r1=octeti(r,1)
      * r2=octeti(r,2)
      * c[0]=sv0^p[0]^r0
      * c[1]=sv1^p[1]^r1
      * c[2]=sv2^p[2]^r2
      */
     int wh=imglin->height*imglin->width;
     uint32_t sv0,sv1,sv2,r0,r1,r2;
     sv0=octeti(sv,0);
     sv1=octeti(sv,1);
     sv2=octeti(sv,2);
    r0=octeti(secv[wh],0);
    r1=octeti(secv[wh],1);
    r2=octeti(secv[wh],2);
    imgtemp2->pixels[0]=(unsigned char)(sv0^imgtemp->pixels[0]^r2);
    imgtemp2->pixels[1]=(unsigned char)(sv1^imgtemp->pixels[1]^r1);
    imgtemp2->pixels[2]=(unsigned char)(sv2^imgtemp->pixels[2]^r0);
    //printf("%d %d %d ",imgtemp2->pixels[0],imgtemp2->pixels[1],imgtemp2->pixels[2]);

    unsigned char cr,cg,cb,pr,pg,pb;
    //restul de pixeli
    for(i=3;i<wh*3;i+=3){
      /*
        c[k]= c[k-1] ^ p'[k] ^ r[wh+k]

        cr[k-1] = c[k-1]
        cg[k-1]= c[k]
        cb[k-1]= c[k+1]

        p'r[k]=p'[k]
        p'g[k]=p'[k+1]
        p'b[k]=p'[k+2]

        r0=octet(r[wh+k],0)
        r1=octet(r[wh+k],1)
        r2=octet(r[wh+k],2)
        */
       //rosu
       cr=imgtemp2->pixels[i-1];
       pr=imgtemp->pixels[i];
       r2=octeti(secv[wh+i],2);
       imgtemp2->pixels[i]=(unsigned char)(cr^pr^r2);
        //verde
        cg=imgtemp2->pixels[i];
        pg=imgtemp->pixels[i+1];
        r1=octeti(secv[wh+i],1);
        imgtemp2->pixels[i+1]=(unsigned char)(cg^pg^r1);
        //rosu
        cb=imgtemp2->pixels[i+1];
        pb=imgtemp->pixels[i+2];
        r0=octeti(secv[wh+i],0);
        imgtemp2->pixels[i+2]=(unsigned char)(cb^pb^r0);

    }
    for(i=0;i<imglin->paddingsize*imglin->height;i++){
        imgtemp2->padding[i]=imglin->padding[i];
    }

    salvare(cript,imgtemp2);
    fclose(f);
    fclose(g);
    fclose(key);
    elibmemorie(imglin);
   elibmemorie(imgtemp);
    elibmemorie(imgtemp2);
    free(state);
    free(secv);
    free(permut);

}

void decriptare(char *init,char *cript,char *cheie){
    FILE *f=fopen(init,"wb");
    FILE *g=fopen(cript,"rb");
    FILE *key=fopen(cheie,"r");
    // incarca imaginea
    struct imgbmp *imglin=incarca(cript);
//nr aleatoare cu xorshift
    uint32_t *state=(uint32_t*)malloc(sizeof(uint32_t));
    uint32_t aux,sv=0;
    fscanf(key,"%d",&aux);
    //printf("%d ",aux);
    fscanf(key,"%d",&sv);
    //printf("%d",sv);

    //fscanf(key,"%d",&sv);
    //printf("%d",sv);
    state[0]=aux;
    uint32_t *secv=(uint32_t*)malloc(((2*imglin->width*imglin->height))*sizeof(uint32_t));
    int i,j;
    for (i = 0; i < (2*imglin->width*imglin->height)-1; i++)
    {
        uint32_t x=xorshift32(state);
        secv[i]=x;
    }
//determinarea permutarii
    uint32_t *permut=(uint32_t*)malloc(imglin->width*imglin->height*sizeof(uint32_t));
    uint32_t p;
    for(p=0;p<imglin->width*imglin->height;p++)
        permut[p]=p;
    for(i=(imglin->height*imglin->width)-1;i>=0;i--){
        j=secv[i]%(i+1);
        aux=permut[i];
        permut[i]=permut[j];
        permut[j]=aux;
    }
//permutarea pixelilor
    struct imgbmp *imgtemp=(struct imgbmp*)malloc(sizeof(struct imgbmp));
    imgtemp->infoheader=(unsigned char*)malloc(54);
    //strcpy(imgtemp->infoheader,imglin->infoheader);
    for(i=0;i<54;i++){
        imgtemp->infoheader[i]=imglin->infoheader[i];
    }
    imgtemp->width = *(unsigned int *)&(imglin->infoheader)[18];
    imgtemp->height = *(unsigned int *)&(imglin->infoheader)[22];
    imgtemp->size=imglin->height*imglin->width*3;
    imgtemp->pixels = malloc(imglin->size);
    imgtemp->paddingsize=imglin->paddingsize;
    imgtemp->padding=(unsigned char*)malloc(imgtemp->paddingsize*imgtemp->height);

    // criptarea in sine
    struct imgbmp *imgtemp2=(struct imgbmp*)malloc(sizeof(struct imgbmp));
    imgtemp2->infoheader=(unsigned char*)malloc(54);
    //strcpy(imgtemp2->infoheader,imglin->infoheader);
    for(i=0;i<54;i++){
        imgtemp2->infoheader[i]=imglin->infoheader[i];
    }
    imgtemp2->width = *(unsigned int *)&(imglin->infoheader)[18];
    imgtemp2->height = *(unsigned int *)&(imglin->infoheader)[22];
    imgtemp2->size=imglin->size;
    imgtemp2->pixels = malloc(imglin->size);
    imgtemp2->paddingsize=imglin->paddingsize;
    imgtemp2->padding=(unsigned char*)malloc(imgtemp2->paddingsize*imgtemp2->height);


    //primul pixel
    /*
     * sv0=octeti(sv,0)
     * sv1=octeti(sv,1)
     * sv2=octeti(sv,2)
     * r0=octeti(r,0)
     * r1=octeti(r,1)
     * r2=octeti(r,2)
     * c[0]=sv0^p[0]^r0
     * c[1]=sv1^p[1]^r1
     * c[2]=sv2^p[2]^r2
     */
    int wh=imglin->height*imglin->width;
    uint32_t sv0,sv1,sv2,r0,r1,r2;
    sv0=octeti(sv,0);
    sv1=octeti(sv,1);
    sv2=octeti(sv,2);
    r0=octeti(secv[wh],0);
    r1=octeti(secv[wh],1);
    r2=octeti(secv[wh],2);
    imgtemp->pixels[0]=(unsigned char)(sv0^imglin->pixels[0]^r2);
    imgtemp->pixels[1]=(unsigned char)(sv1^imglin->pixels[1]^r1);
    imgtemp->pixels[2]=(unsigned char)(sv2^imglin->pixels[2]^r0);
    //printf("%d %d %d ",imgtemp2->pixels[0],imgtemp2->pixels[1],imgtemp2->pixels[2]);

    unsigned char cr,cg,cb,pr,pg,pb;
    //restul de pixeli
    for(i=3;i<wh*3;i+=3){
        /*
         c[k]= c[k-1] ^ p'[k] ^ r[wh+k]

         cr[k-1] = c[k-1]
         cg[k-1]= c[k]
         cb[k-1]= c[k+1]

         p'r[k]=p'[k]
         p'g[k]=p'[k+1]
         p'b[k]=p'[k+2]

         r0=octet(r[wh+k],0)
         r1=octet(r[wh+k],1)
         r2=octet(r[wh+k],2)
         */
        //rosu
        cr=imglin->pixels[i-1];
        pr=imglin->pixels[i];
        r2=octeti(secv[wh+i],2);
        imgtemp->pixels[i]=(unsigned char)(cr^pr^r2);
        //verde
        cg=imglin->pixels[i];
        pg=imglin->pixels[i+1];
        r1=octeti(secv[wh+i],1);
        imgtemp->pixels[i+1]=(unsigned char)(cg^pg^r1);
        //rosu
        cb=imglin->pixels[i+1];
        pb=imglin->pixels[i+2];
        r0=octeti(secv[wh+i],0);
        imgtemp->pixels[i+2]=(unsigned char)(cb^pb^r0);

    }

    for(i=0;i<imgtemp->width*imgtemp->height;i++){
        {imgtemp2->pixels[i*3]=imgtemp->pixels[permut[i]*3];
            imgtemp2->pixels[i*3+1]=imgtemp->pixels[permut[i]*3+1];
            imgtemp2->pixels[i*3+2]=imgtemp->pixels[permut[i]*3+2];
        }
    }

   for(i=0;i<imglin->paddingsize*imglin->height;i++)
       imgtemp2->padding[i]=imglin->padding[i];

    salvare(init,imgtemp2);
    elibmemorie(imglin);
    elibmemorie(imgtemp);
    elibmemorie(imgtemp2);
    free(state);
    free(secv);
    free(permut);
    fclose(f);
    fclose(g);
    fclose(key);
}

void testchi(char *str) {
        int x,i,j;
        float rsum=0,gsum=0,bsum=0,med;
        struct imgbmp *imglin=incarca(str);
        med=((float)(imglin->width*imglin->height)/256);
        printf("\nmed=%f\n",med);
        float *fr=(float*)malloc(256 * sizeof(float));
        float *fg=(float*)malloc(256 * sizeof(float));
        float *fb=(float*)malloc(256 * sizeof(float));
        for(i=0;i<256;i++)
            fr[i]=fg[i]=fb[i]=0;
        for(i=0;i<imglin->width*imglin->height;i++){
            fr[imglin->pixels[i*3]]++;
            fg[imglin->pixels[i*3+1]]++;
            fb[imglin->pixels[i*3+2]]++;
        }
        for(i=0;i<256;i++){
            rsum+=((fr[i]-med)*(fr[i]-med))/med;
            gsum+=((fg[i]-med)*(fg[i]-med))/med;
            bsum+=((fb[i]-med)*(fb[i]-med))/med;
        }
        printf("\nValoarea pentru canalul R este: %f",rsum);
        printf("\nValoarea pentru canalul G este: %f",gsum);
        printf("\nValoarea pentru canalul B este: %f",bsum);
        elibmemorie(imglin);
        free(fr);
        free(fg);
        free(fb);
}


void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], header[54], aux;


    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL)
    {
        printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    //copiaza octet cu octet imaginea initiala in cea noua
    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);

    //calculam padding-ul pentru o linie
    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;


    fseek(fout, 54, SEEK_SET);
    int i,j;
    for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            //citesc culorile pixelului
            fread(pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pRGB, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fout);
}

struct coord* tempmatch(char *img,char *sablon,float ps){
    //FILE *f=fopen(img,"rb");
    //FILE *g=fopen(sablon,"rb");
    //grayscale_image(img,"imggray.bmp");
    //grayscale_image(sablon,"sablongray.bmp");
    struct imgbmp* imggray=incarca(img);
    struct imgbmp* sabgray=incarca(sablon);
    //struct coord2 *ferestre2=(struct coord2*)malloc(sizeof(struct coord2)*imggray->height*imggray->width);

    struct coord *ferestre=(struct coord*)malloc(sizeof(struct coord));
    ferestre->x=(long*)malloc(sizeof(long)*imggray->height*imggray->width);
    ferestre->y=(long*)malloc(sizeof(long)*imggray->height*imggray->width);
    ferestre->corr=(float*)malloc(sizeof(float)*imggray->height*imggray->width);
    ferestre->lungime=0;
    int i,j,k,l;
    float smediu=0,suma=0;
    float fImediu=0,sigmas=0,sigmafI=0,corr=0;
    //suma medie sablon
for(i=0;i<sabgray->height;i++)
  {
            for(j=0;j<sabgray->width;j++)
                {
                    suma+=(float)sabgray->pixels[i*sabgray->width*3+j*3];
                }
     }

     smediu=suma/(sabgray->height*sabgray->width);
//printf("%f",smediu);
//sigma sablon
    for(i=0;i<sabgray->height;i++)
            {
                for(j=0;j<sabgray->width;j++)
                {
                    sigmas+=(((float)(sabgray->pixels[i*sabgray->width*3+j*3]))-smediu)*(((float)(sabgray->pixels[i*sabgray->width*3+j*3]))-smediu);

                }
            }
        sigmas=sigmas/(sabgray->height*sabgray->width-1);
        //printf("\n%f",sigmas);
        sigmas=sqrt(sigmas);
        //printf("\n%f",sigmas);
        //fImediu , sigma fereastra , corelatia
	for(i=0;i<imggray->height-sabgray->height+1;i++){
        for(j=0;j<imggray->width-sabgray->width+1;j++)
        {
            sigmafI=0;
            fImediu=0;
            for(k=0;k<sabgray->height;k++)
            {
                for(l=0;l<sabgray->width;l++)
                {
                    fImediu+=imggray->pixels[i*imggray->width*3+j*3+k*imggray->width*3+l*3];
                }
            }
            //printf("\n%d",i);
            fImediu=fImediu/(sabgray->height*sabgray->width);
            for(k=0;k<sabgray->height;k++)
            {
                for(l=0;l<sabgray->width;l++)
                {
                    sigmafI+=(imggray->pixels[i*imggray->width*3+j*3+k*imggray->width*3+l*3]-fImediu)*(imggray->pixels[i*imggray->width*3+j*3+k*imggray->width*3+l*3]-fImediu);
                }
            }
            for(k=0;k<sabgray->height;k++)
            {
                for(l=0;l<sabgray->width;l++)
                {
                    corr+=(imggray->pixels[i*imggray->width*3+j*3+k*imggray->width*3+l*3]-fImediu)*(sabgray->pixels[k*sabgray->width*3+l*3]-smediu);
                }

            }

            sigmafI=sigmafI/(sabgray->width*sabgray->height-1);
    sigmafI=sqrt(sigmafI);
    corr=((corr/sigmafI)/sigmas)/(sabgray->width*sabgray->height);
    //printf("%f,",sigmas);
            //printf("corr=%f",corr);
            //printf("\n");
    //daca corelatia>ps,se salveaza coordonatele si valoarea corelatiei

    if(corr>ps)
    {
    ferestre->corr[ferestre->lungime]=corr;
    ferestre->x[ferestre->lungime]=i;
    ferestre->y[ferestre->lungime]=j;
    ferestre->lungime++;
        //printf("corr=%f",corr);
    //printf("\n");
     }
    corr=0;
    fImediu=0;
        }
	}
	//printf("\n%f %f",smed,fImediu);
	//printf("%d",cont);
    ferestre->x=realloc(ferestre->x,sizeof(long)*ferestre->lungime);
    ferestre->y=realloc(ferestre->y,sizeof(long)*ferestre->lungime);
    ferestre->corr=realloc(ferestre->corr,sizeof(float)*ferestre->lungime);

    elibmemorie(imggray);
    elibmemorie(sabgray);
    //fclose(f);
    //fclose(g);
    return ferestre;
}

void coloreazafereastra(char *img,char *sablonim, long x, long y,int cifra){
struct imgbmp *imglin=incarca(img);
struct imgbmp *sablon=incarca(sablonim);
int i,j;
unsigned char cr,cg,cb;
//culori pentru fiecare cifra
    if(cifra==0){
        cr=255;
        cg=0;
        cb=0;
    }
    if(cifra==1){
    cr=255;
    cg=255;
    cb=0;
    }
    if(cifra==2){
        cr=0;
        cg=255;
        cb=0;
    }
    if(cifra==3){
        cr=0;
        cg=255;
        cb=255;
    }
    if(cifra==4){
        cr=255;
        cg=0;
        cb=255;
    }
    if(cifra==5){
        cr=0;
        cg=0;
        cb=255;
    }
    if(cifra==6){
        cr=192;
        cg=192;
        cb=192;
    }
    if(cifra==7){
        cr=255;
        cg=140;
        cb=0;
    }
    if(cifra==8){
        cr=128;
        cg=0;
        cb=128;
    }
    if(cifra==9){
        cr=128;
        cg=0;
        cb=0;
    }

//m am plictisit de proiect
//colorarea imaginii

if(x<=imglin->height-sablon->height && y<=imglin->width-sablon->width) {
   for (i = 0; i < sablon->width; i++) {
        imglin->pixels[(x * imglin->width + y + i) * 3] = cb;
        imglin->pixels[(x * imglin->width + y + i) * 3 + 1] = cg;
        imglin->pixels[(x * imglin->width + y + i) * 3 + 2] = cr;
    }

    for (i = 0; i <= sablon->width; i++) {
        imglin->pixels[((x + sablon->height) * imglin->width + y + i) * 3] = cb;
        imglin->pixels[((x + sablon->height) * imglin->width + y + i) * 3 + 1] = cg;
        imglin->pixels[((x + sablon->height) * imglin->width + y + i) * 3 + 2] = cr;
    }

    for (i = 0; i < sablon->height; i++) {
        imglin->pixels[((x + i) * imglin->width + y) * 3] = cb;
        imglin->pixels[((x + i) * imglin->width + y) * 3 + 1] = cg;
        imglin->pixels[((x + i) * imglin->width + y) * 3 + 2] = cr;
    }

    for (i = 0; i < sablon->height; i++) {
        imglin->pixels[((x + i) * imglin->width + y + sablon->width) * 3] = cb;
        imglin->pixels[((x + i) * imglin->width + y + sablon->width) * 3 + 1] = cg;
        imglin->pixels[((x + i) * imglin->width + y + sablon->width) * 3 + 2] = cr;
    }

}

salvare(img,imglin);
elibmemorie(imglin);
elibmemorie(sablon);
}

int cmp(const void *a, const void *b){
    struct coord2 va=*(struct coord2*)a;
    struct coord2 vb=*(struct coord2*)b;
    if(va.corr<vb.corr)
        return 1;
    else
        return -1;
}

void sortarecorelatii(struct coord **corr){
    int i;
    struct coord2 *ferestre2=(struct coord2*)malloc(sizeof(struct coord2)*(*corr)->lungime);
    for(i=0;i<(*corr)->lungime;i++){
        ferestre2[i].x=(*corr)->x[i];
        ferestre2[i].y=(*corr)->y[i];
        ferestre2[i].corr=(*corr)->corr[i];
        ferestre2[i].cifra=(*corr)->cifra[i];
    }

    qsort(ferestre2,(*corr)->lungime,sizeof(struct coord2),cmp);

    for(i=0;i<(*corr)->lungime;i++){
        (*corr)->x[i]=ferestre2[i].x;
        (*corr)->y[i]=ferestre2[i].y;
        (*corr)->corr[i]=ferestre2[i].corr;
        (*corr)->cifra[i]=ferestre2[i].cifra;
    }

    free(ferestre2);
}


void eliminarenonmaxime(struct coord **corr,char *str){
    long i,j,x,y,v,w,k;
    float aij;
    float suprapunere;
    struct imgbmp* sablon=incarca(str);
    float arie=sablon->width*sablon->height;
    for(i=0;i<(*corr)->lungime-1;i++){
        for(j=i+1;j<(*corr)->lungime;j++){
            if(labs((*corr)->x[i]-(*corr)->x[j])<=sablon->height && labs((*corr)->y[i]-(*corr)->y[j])<=sablon->width){
                w=labs((*corr)->x[i]-(*corr)->x[j]);
                x=labs(sablon->height-w);
                v=labs((*corr)->y[i]-(*corr)->y[j]);
                y=labs(sablon->width-v);
                aij=x*y;
                suprapunere=(aij/(arie+arie+aij));
                if(suprapunere>=0.2){
                    for(k=j;k<(*corr)->lungime-1;k++){
                        (*corr)->x[k]=(*corr)->x[k+1];
                        (*corr)->y[k]=(*corr)->y[k+1];
                        (*corr)->corr[k]=(*corr)->corr[k+1];
                        (*corr)->cifra[k]=(*corr)->cifra[k+1];
                    }
                    (*corr)->lungime--;
                    j=0;
                    i=0;
                }
            }
        }
    }

    (*corr)->x=realloc((*corr)->x,(*corr)->lungime*sizeof(long));
    (*corr)->y=realloc((*corr)->y,(*corr)->lungime*sizeof(long));
    (*corr)->corr=realloc((*corr)->corr,(*corr)->lungime*sizeof(float));
    (*corr)->cifra=realloc((*corr)->cifra,(*corr)->lungime* sizeof(int));
    elibmemorie(sablon);

}



int main()
{
//criptare,decriptare

    char *pathimg=(char*)malloc(200);
    char *pathsablon=(char*)malloc(200);
    char *pathcript=(char*)malloc(200);
    char *pathkey=(char*)malloc(200);
    char *pathimginit=(char*)malloc(200);
    printf("\nCalea imaginii de criptat:\n");
    scanf("%s",pathimginit);
    printf("\nCalea unde va fi salvata imaginea criptata:\n");
    scanf("%s",pathcript);
    printf("\nCalea cheii criptarii");
    scanf("%s",pathkey);
        criptare(pathimginit,pathcript,pathkey);
    printf("\nCalea imaginii de decriptat:\n");
    scanf("%s",pathcript);
    printf("\nCalea unde va fi salvata imaginea decriptata:\n");
    scanf("%s",pathimg);
    printf("\nCalea cheii decriptarii:\n");
    scanf("%s",pathkey);
        decriptare(pathimg,pathcript,pathkey);
    printf("\nTestul chi pentru imaginea initiala:\n")
        testchi(pathimginit);
    printf("\nTestul chi pentru imaginea criptata:\n")
        testchi(pathcript);


//recunoasterea cifrelor


    int nrsabloane,j,i;
    float ps=0.5;
    struct coord *corrtotal=(struct coord*)malloc(sizeof(struct coord));
   printf("\nCalea imaginii pentru template matching:\n");
   scanf("%s",pathimg);
   grayscale_image(pathimg,"imaginegray.bmp");
   printf("\nNumar sabloane:\n");
   scanf("%d",&nrsabloane);
   for(i=0;i<nrsabloane;i++){
    printf("\nCalea sablonului numarul %d:\n",i);
    scanf("%s",pathsablon);
    grayscale_image(pathsablon,"sablongray.bmp");
    struct coord *corrpartial=tempmatch("imaginegray.bmp","sablongray.bmp",ps);
    if(i==0){
        corrtotal->x=(long*)malloc(sizeof(long)*corrpartial->lungime);
        corrtotal->y=(long*)malloc(sizeof(long)*corrpartial->lungime);
        corrtotal->corr=(float*)malloc(sizeof(float)*corrpartial->lungime);
        corrtotal->cifra=(int*)malloc(sizeof(int)*corrpartial->lungime);
        corrtotal->lungime=corrpartial->lungime;
        for(j=0;j<corrpartial->lungime;j++){
            corrtotal->x[j]=corrpartial->x[j];
            corrtotal->y[j]=corrpartial->y[j];
            corrtotal->corr[j]=corrpartial->y[j];
            corrtotal->cifra[j]=0;
        }
    }
    else {

        corrtotal->x=realloc(corrtotal->x,sizeof(long)*(corrtotal->lungime+corrpartial->lungime));
        corrtotal->y=realloc(corrtotal->y,sizeof(long)*(corrtotal->lungime+corrpartial->lungime));
        corrtotal->corr=realloc(corrtotal->corr,sizeof(float)*(corrtotal->lungime+corrpartial->lungime));
        corrtotal->cifra=realloc(corrtotal->cifra,sizeof(int)*(corrtotal->lungime+corrpartial->lungime));
        for(j=0;j<corrpartial->lungime;j++){
            corrtotal->x[j+corrtotal->lungime]=corrpartial->x[j];
            corrtotal->y[j+corrtotal->lungime]=corrpartial->y[j];
            corrtotal->corr[j+corrtotal->lungime]=corrpartial->corr[j];
            corrtotal->cifra[j+corrtotal->lungime]=i;
        }
        corrtotal->lungime+=corrpartial->lungime;
    }

        remove("sablongray.bmp");
   }

   sortarecorelatii(&corrtotal);
   eliminarenonmaxime(&corrtotal,pathsablon);
   for(i=0;i<corrtotal->lungime;i++){
    coloreazafereastra(pathimg,pathsablon,corrtotal->x[i],corrtotal->y[i],corrtotal->cifra[i]);
   }
   elibmemorie2(corrtotal);
   free(pathimg);
   free(pathsablon);
   free(pathcript);
   free(pathkey);
   free(pathimginit);
}