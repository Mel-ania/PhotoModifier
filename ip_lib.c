/*
 * Bonesso Max
 * Gottardo Melania
 */

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l=0; l<t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0; i<t->h; i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}


/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned  int k, float v){
  unsigned int height=0, width=0, channel=0;
  ip_mat *a = (ip_mat*)malloc(sizeof(ip_mat));
  if(!a){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(a);
    exit(1);
  }
  a->h = (float) h;
  a->w = (float) w;
  a->k = (float) k;
  a->stat = (stats*)malloc(sizeof(stats)*k);
  if(!a->stat){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(a);
    exit(1);
  }
  for(channel=0; channel<k; channel++){
    a->stat[channel].min = v;
    a->stat[channel].max = v;
    a->stat[channel].mean = v;
  }
  a->data = (float***)malloc(sizeof(float**)*h);
  if(!a->data){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(a);
    exit(1);
  }
  for(height=0; height<h; height++){
    a->data[height] = (float**)malloc(sizeof(float*)*w);
    if(!a->data[height]){
      printf("Errore: allocazione memoria\n");
      ip_mat_free(a);
      exit(1);
    }
  }
  for(height=0; height<h; height++){
    for(width=0; width<w; width++){
      a->data[height][width] = (float*)malloc(sizeof(float)*k);
      if(!a->data[height][width]){
        printf("Errore: allocazione memoria\n");
        ip_mat_free(a);
        exit(1);
      }
      for(channel=0; channel<k; channel++){
        a->data[height][width][channel] = v;
      }
    }
  }
  return a;
}

void ip_mat_free(ip_mat *a){
  unsigned int height=0, width=0;
  if(a){
    if(a->data){
      for(height=0; height<a->h; height++){
        if(a->data[height]){
          for(width=0; width<a->w; width++){
            if(a->data[height][width])
              free(a->data[height][width]);
            }
          free(a->data[height]);
        }
      }
      free(a->data);
    }
    if(a->stat)
      free(a->stat);
    free(a);
  }
}

void compute_stats(ip_mat * t){
  unsigned int channel=0, width=0, height=0;
  int count=0;
  float num=0., res=0.;
  for(channel=0; channel<t->k; channel++){
    t->stat[channel].min = t->data[0][0][channel];
    t->stat[channel].max = t->data[0][0][channel];
  }
  for(channel=0; channel<t->k; channel++){
	for(height=0; height<t->h; height++){
		for(width=0; width<t->w; width++){
			num = t->data[height][width][channel];
			if(t->stat[channel].min > num){
			  t->stat[channel].min = num;
			}
			if(t->stat[channel].max < num){
			  t->stat[channel].max = num;
			}
			count++;
			res+=num;
      }
    }
    t->stat[channel].mean = res/(float)count;
    count=0;
    res=0;
  }
}

void ip_mat_init_random(ip_mat * t, float mean, float var){
  unsigned int height=0, width=0, channel=0;
  for(height=0; channel<t->h; height++){
    for(width=0; width<t->w; width++){
      for(channel=0; channel<t->k; channel++){
        t->data[height][width][channel] = get_normal_random(mean, var);
      }
    }
  }
  compute_stats(t);
}

ip_mat * ip_mat_copy(ip_mat * in){
  unsigned int channel=0, width=0, height=0;
  ip_mat *out = (ip_mat*)malloc(sizeof(ip_mat));
  if(!out){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  out->h = in->h;
  out->w = in->w;
  out->k = in->k;
  out->stat = (stats*)malloc(sizeof(stats)*(in->k));
  if(!out->stat){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  for(channel=0; channel<in->k; channel++){
    out->stat[channel].min = in->stat[channel].min;
    out->stat[channel].max = in->stat[channel].max;
    out->stat[channel].mean = in->stat[channel].mean;
  }
  out->data = (float***)malloc(sizeof(float**)*(in->h));
  if(!out->data){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  for(height=0; height<in->h; height++){
    out->data[height] = (float**)malloc(sizeof(float*)*(in->w));
    if(!out->data[height]){
      printf("Errore: allocazione memoria\n");
      ip_mat_free(out);
      exit(1);
    }
  }
  for(height=0; height<in->h; height++){
    for(width=0; width<in->w; width++){
      out->data[height][width] = (float*)malloc(sizeof(float)*(in->k));
      if(!out->data[height][width]){
        printf("Errore: allocazione memoria\n");
        ip_mat_free(out);
        exit(1);
      }
      for(channel=0; channel<in->k; channel++){
        out->data[height][width][channel] = in->data[height][width][channel];
      }
    }
  }
  return out;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
  unsigned int channel=0, width=0, height=0;
  ip_mat *out = NULL;
  if((t->h)<row_end-1 || row_start>=row_end){
    printf("Errore: overflow dimensioni righe\n");
    ip_mat_free(out);
    exit(1);
  }
  if((t->w)<col_end-1 || col_start>=col_end){
    printf("Errore: overflow dimensioni colonne\n");
    ip_mat_free(out);
    exit(1);
  }
  out = (ip_mat*)malloc(sizeof(ip_mat));
  if(!out){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  out->w = col_end-col_start;
  out->h = row_end-row_start;
  out->k = t->k;
  out->stat = (stats*)malloc(sizeof(stats)*(out->k));
  if(!out->stat){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  out->data = (float***)malloc(sizeof(float**)*(out->h));
  if(!out->data){
    printf("Errore: allocazione memoria\n");
    ip_mat_free(out);
    exit(1);
  }
  for(height=0; height<out->h; height++){
    out->data[height] = (float**)malloc(sizeof(float*)*(out->w));
    if(!out->data[height]){
      printf("Errore: allocazione memoria\n");
      ip_mat_free(out);
      exit(1);
    }
  }
  for(height=0; height<out->h; height++){
    for(width=0; width<out->w; width++){
      out->data[height][width] = (float*)malloc(sizeof(float)*(out->k));
      if(!out->data[height][width]){
        printf("Errore: allocazione memoria\n");
        ip_mat_free(out);
        exit(1);
      }
      for(channel=0; channel<out->k; channel++){
        out->data[height][width][channel] = t->data[height+row_start][width+col_start][channel];
      }
    }
  }
  compute_stats(out);
  return out;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
  ip_mat *out = NULL;
  unsigned int height=0, width=0, channel=0;
  if(dimensione == 0){
    if(a->w!=b->w || a->k!=b->k) {
      printf("Errore: dimensioni colonne e canali diverse\n");
      ip_mat_free(out);
      exit(1);
    }
    out = ip_mat_create(a->h+b->h, a->w, a->k, 0.);
    for(height=0; height<a->h+b->h; height++){
      for(width=0; width<a->w; width++){
        for(channel=0; channel<a->k; channel++){
          out->data[height][width][channel] = (height < a->h) ? (a->data[height][width][channel]) : (b->data[height-a->h][width][channel]);
        }
      }
    }
  }

  else if(dimensione == 1){
    if(a->h!=b->h || a->k!=b->k){
      printf("Errore: dimensioni righe e canali diverse\n");
      ip_mat_free(out);
      exit(1);
    }
    out = ip_mat_create(a->h, a->w+b->w, a->k, 0.);
    for(height=0; height<a->h; height++){
      for(width=0; width<a->w+b->w; width++){
        for(channel=0; channel<a->k; channel++){
          out->data[height][width][channel] = (width < a->w) ? (a->data[height][width][channel]) : (b->data[height][width-a->w][channel]);
        }
      }
    }
  }

  else if(dimensione == 2){
    if(a->h!=b->h || a->w!=b->w){
      printf("Errore: dimensioni colonne e righe diverse\n");
      ip_mat_free(out);
      exit(1);
    }
    out = ip_mat_create(a->h, a->w, a->k+b->k, 0.);
    for(height=0; height<a->h; height++){
      for(width=0; width<a->w; width++){
        for(channel=0; channel<a->k+b->k; channel++){
          out->data[height][width][channel] = (channel < a->k) ? (a->data[height][width][channel]) : (b->data[height][width][channel-a->k]);
        }
      }
    }
  }

  else {
    printf("Errore: dimensione non supportata\n");
    ip_mat_free(out);
    exit(1);
  }
  compute_stats(out);
  return out;
}


/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
  ip_mat * sum = NULL;
  unsigned int height=0, width=0, channel=0;
  if(a->h != b->h || a->w != b->w || a->k != b->k){
    printf("Errore: dimensioni diverse");
    ip_mat_free(sum);
    exit(1);
  }
  sum = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        sum->data[height][width][channel] += b->data[height][width][channel];
      }
    }
  }
  compute_stats(sum);
  return sum;
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
  ip_mat * sub = NULL;
  unsigned int height=0, width=0, channel=0;
  if(a->h != b->h || a->w != b->w || a->k != b->k){
    printf("Errore: dimensioni diverse");
    ip_mat_free(sub);
    exit(1);
  }
  sub = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        sub->data[height][width][channel] -= b->data[height][width][channel];
      }
    }
  }
  compute_stats(sub);
  return sub;
}

ip_mat * ip_mat_mul(ip_mat * a, ip_mat * b){
  ip_mat * mul = NULL;
  unsigned int height=0, width=0, channel=0;
  if(a->h != b->h || a->w != b->w || a->k != b->k){
    printf("Errore: dimensioni diverse");
    ip_mat_free(mul);
    exit(1);
  }
  mul = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        mul->data[height][width][channel] *= b->data[height][width][channel];
      }
    }
  }
  compute_stats(mul);
  return mul;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
  ip_mat * ms = NULL;
  unsigned int height=0, width=0, channel=0;
  ms = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        ms->data[height][width][channel] *= c;
      }
    }
  }
  compute_stats(ms);
  return ms;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
  ip_mat * as = NULL;
  unsigned int height=0, width=0, channel=0;
  as = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        as->data[height][width][channel] += c;
      }
    }
  }
  compute_stats(as);
  return as;
}

float ip_mat_sum_mat(ip_mat * a, int channel){
  float res=0.;
  unsigned int height=0, width=0;
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      res+=a->data[height][width][channel];
    }
  }
  return res;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
  ip_mat * mn = NULL;
  unsigned int height=0, width=0, channel=0;
  if(a->h != b->h || a->w != b->w || a->k != b->k){
    printf("Errore: dimensioni diverse");
    ip_mat_free(mn);
    exit(1);
  }
  mn = ip_mat_copy(a);
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        mn->data[height][width][channel] += b->data[height][width][channel];
        mn->data[height][width][channel] /= 2;
      }
    }
  }
  compute_stats(mn);
  return mn;
}


/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
  ip_mat * gray = NULL;
  unsigned int height=0, width=0, channel=0;
  float tot=0.;
  gray = ip_mat_copy(in);
  for(height=0; height<in->h; height++){
    for(width=0; width<in->w; width++){
      for(channel=0; channel<in->k; channel++){
        tot += gray->data[height][width][channel];
      }
      for(channel=0; channel<in->k; channel++){
        gray->data[height][width][channel] = tot/in->k;
      }
      tot = 0.;
    }
  }
  compute_stats(gray);
  return gray;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
  ip_mat * blend = NULL;
  unsigned int height=0, width=0, channel=0;
  blend = ip_mat_create(a->h, a->w, a->h, 0);
  if(a->h!=b->h || a->w!=b->w || a->k!=b->k){
    printf("Errore: dimensioni diverse");
    ip_mat_free(blend);
    exit(1);
  }
  for(height=0; height<a->h; height++){
	  for(width=0; width<a->w; width++){
		  for(channel=0; channel<a->k; channel++){
			  blend->data[height][width][channel] = alpha * a->data[height][width][channel] + (1-alpha) * b->data[height][width][channel];
          }
      }
  }
  compute_stats(blend);
  return blend;
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
  return ip_mat_add_scalar(a, bright);
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
  ip_mat * crpt = ip_mat_copy(a);
  unsigned int height=0, width=0, channel=0;
  for(height=0; height<a->h; height++){
    for(width=0; width<a->w; width++){
      for(channel=0; channel<a->k; channel++){
        crpt->data[height][width][channel] = a->data[height][width][channel] + amount*get_normal_random(0,1);
      }
    }
  }
  compute_stats(crpt);
  return crpt;
}


/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
  ip_mat * n = ip_mat_padding(a, (f->h-1)/2, (f->w-1)/2);
  ip_mat * c = ip_mat_copy(a);
  ip_mat * sup = NULL, *sup2 = NULL, *f2 = ip_mat_copy(f);
  unsigned int height=0, width=0, channel=0;
  while(a->k > f2->k){
    ip_mat * f1 = ip_mat_copy(f2);
    ip_mat_free(f2);
    f2 = ip_mat_concat(f1, f, 2);
    ip_mat_free(f1);
  }
  for(height=0; height<(n->h)-(f2->h); height++){
    for(width=0; width<(n->w)-(f2->w); width++){
		sup = ip_mat_subset(n, height, height + f2->h, width, width + f2->w);
		sup2 = ip_mat_mul(sup, f2);
		for(channel=0; channel<n->k; channel++){
			c->data[height][width][channel] = ip_mat_sum_mat(sup2, channel);
		}
        ip_mat_free(sup);
        ip_mat_free(sup2);
    }
  }
  compute_stats(c);
  ip_mat_free(n);
  ip_mat_free(f2);
  return c;
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
  unsigned int height=0, width=0, channel=0;
  ip_mat * n;
  n = ip_mat_create((a->h)+2*pad_h, (a->w)+2*pad_w, a->k, 0);
  for(height=pad_h; height<(a->h)+pad_h-1; height++){
    for(width=pad_w; width<(a->w)+pad_w-1; width++){
      for(channel=0; channel<a->k; channel++){
        n->data[height][width][channel] = a->data[height-pad_h][width-pad_w][channel];
      }
    }
  }
  compute_stats(n);
  return n;
}

ip_mat * create_sharpen_filter(){
  ip_mat * s = ip_mat_create(3, 3, 1, 0);
  s->data[0][1][0]=-1;
  s->data[1][0][0]=-1;
  s->data[1][1][0]=5;
  s->data[1][2][0]=-1;
  s->data[2][1][0]=-1;
  compute_stats(s);
  return s;
}

ip_mat * create_edge_filter(){
  ip_mat * e = ip_mat_create(3, 3, 1, -1);
  e->data[1][1][0]=8;
  compute_stats(e);
  return e;
}

ip_mat * create_emboss_filter(){
  ip_mat * em = ip_mat_create(3, 3, 1, 1);
  em->data[0][0][0]=-2;
  em->data[0][1][0]=-1;
  em->data[0][2][0]=0;
  em->data[1][0][0]=-1;
  em->data[2][0][0]=0;
  em->data[2][2][0]=2;
  compute_stats(em);
  return em;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
  return ip_mat_create(h, w, k, (1./(float)(w*h)));
}

ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma){
    ip_mat * g = ip_mat_create(h, w, k, 0);
    int cx = h/2;
    int cy = w/2;
    unsigned int width=0, height=0, channel=0;
    int x=0, y=0;
    float sum = 0.;
    for(height=0; height<h; height++){
      for(width=0; width<w; width++){
        for(channel=0; channel<k; channel++){
          y = width - cy;
          x = height - cx;
          g->data[height][width][channel] = (1/(2*PI*sigma*sigma))*exp(-1*((x*x+y*y)/(2*sigma*sigma)));
        }
      }
    }
    for(channel=0; channel<k; channel++){
      ip_mat * g1 = ip_mat_copy(g);
      sum = ip_mat_sum_mat(g, channel);
      ip_mat_free(g);
      g = ip_mat_mul_scalar(g1, 1/sum);
      sum = 0;
      ip_mat_free(g1);
    }
    compute_stats(g);
    return g;
}

void rescale(ip_mat * t, float new_max){
  unsigned int height=0, width=0, channel=0;
  for(height=0; height<t->h; height++){
    for(width=0; width<t->w; width++){
      for(channel=0; channel<t->k; channel++){
        t->data[height][width][channel] = (t->data[height][width][channel] - t->stat[channel].min)/(t->stat[channel].max - t->stat[channel].min);
      }
    }
  }
  t = ip_mat_mul_scalar(t, new_max);
  compute_stats(t);
}

void clamp(ip_mat * t, float low, float high){
  unsigned int height=0, width=0, channel=0;
  for(height=0; height<t->h; height++){
    for(width=0; width<t->w; width++){
      for(channel=0; channel<t->k; channel++){
        if(t->data[height][width][channel] > high){
          t->data[height][width][channel] = high;
        }
        if(t->data[height][width][channel] < low){
          t->data[height][width][channel] = low;
        }
      }
    }
  }
  compute_stats(t);
}
