/*
 * CYANCORE LICENSE
 * Copyrights (C) 2022, Cyancore Team
 *
 * File Name		: malloc_lite.c
 * Description		: This file contains sources of libc-malloc
 *			  functions
 * Primary Author	: Akash Kollipara [akashkollipara@gmail.com]
 * Organisation		: Cyancore Core-Team
 */

#include <stdint.h>
// #include <status.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "nnom_port.h"
#include <stdarg.h>
#include <stdbool.h>

#define get_num_va_args(_args, _lcount)			\
	(((_lcount) >= 2) ? va_arg(_args, int64_t) :	\
	((_lcount) == 1) ? va_arg(_args, long) :	\
			    va_arg(_args, int))

#define get_unum_va_args(_args, _lcount)			\
	(((_lcount) >= 2) ? va_arg(_args, uint64_t) :		\
	((_lcount) == 1) ? va_arg(_args, unsigned long) :	\
			    va_arg(_args, unsigned int))
#define USE_FLOAT 1
// #include <lock/lock.h>
// #include <arch.h>
// #include <mmio.h>
// #include <plat_mem.h>
// #define HEAP_SIZE 16K
#define HEAP_ALIGN 4
#define ROUNDUP_ALIGN(x, align)        x += align - (x % align)

extern uint8_t _heap_start, _heap_end, _heap_size;
// static istate_t state;
// static lock_t mlock;

typedef struct chunk
{
	size_t size;
	size_t free;
	struct chunk *next;
} chunk_t;

static chunk_t *freeList;

// static void heap_lock(void)
// {
// 	lock_acquire(&mlock);
// 	arch_di_save_state(&state);
// }

// static void heap_unlock(void)
// {
// 	lock_release(&mlock);
// 	arch_ei_restore_state(&state);
// }

static void split(chunk_t *fit_slot, size_t size)
{
	chunk_t *new = (chunk_t *)((size_t) fit_slot + size + sizeof(chunk_t));
	new->size = fit_slot->size - size - sizeof(chunk_t);
	new->free = 1;
	new->next = fit_slot->next;
	fit_slot->size = size;
	fit_slot->free = 0;
	fit_slot->next = new;
	return;
}

static void merge()
{
	chunk_t *cur;
	cur = freeList;
	while(cur && cur->next)
	{
		if(cur->free && cur->next->free)
		{
			cur->size += cur->next->size + sizeof(chunk_t);
			cur->next = cur->next->next;
		}
		cur = cur->next;
	}
}

static chunk_t *get_header(void *p)
{
	return (chunk_t *)p - 1;
}

// status_t platform_init_heap()
// {
// 	heap_lock();
// 	size_t sz = (size_t)&_heap_end;
// 	sz -= (size_t)&_heap_start;
// 	memset(&_heap_start, 0, sz);
// 	freeList = (void *)&_heap_start;
// 	freeList->size = sz - sizeof(chunk_t);
// 	freeList->free = 1;
// 	freeList->next = NULL;
// 	heap_unlock();
// 	return success;
// }
void platform_init_heap()
{
	// heap_lock();
	size_t sz = (size_t)&_heap_end;
	sz -= (size_t)&_heap_start;
	memset(&_heap_start, 0, sz);
	freeList = (void *)&_heap_start;
	freeList->size = sz - sizeof(chunk_t);
	freeList->free = 1;
	freeList->next = NULL;
	// heap_unlock();
	return 0;
}
void *my_malloc(size_t n_bytes)
{
	chunk_t *cur;
	void *p = NULL;

	if(!n_bytes)
		return NULL;


	ROUNDUP_ALIGN(n_bytes, 4);

	// heap_lock();
	if(!freeList->size)
		goto exit;
	cur = freeList;
	while(cur->size < n_bytes || (!cur->free && cur->next))
		cur = cur->next;
	if(cur->size == n_bytes)
		cur->free = 0;
	else if(cur->size > (n_bytes + sizeof(chunk_t)))
		split(cur, n_bytes);
	else
		goto exit;
	p = (void *)(++cur);
exit:
	// heap_unlock();
	return p;
}

int  my_free(void *ptr)
{
	if(ptr == NULL)
		return 0;

	chunk_t *cur = get_header(ptr);
	// heap_lock();
	if((void *)&_heap_start <= (void *)cur &&
		(void *)cur <= (void *)&_heap_end)
	{
		cur->free = 1;
		merge();
		merge();
		(volatile void ) 0;
	}
	// heap_unlock();
	return 0;
}
#if 0
void *calloc(size_t n_blocks, size_t n_bytes)
{
	n_bytes *= n_blocks;
	void *p = malloc(n_bytes);
	if(p)
		memset(p, 0, n_bytes);
	return p;
}

void *realloc(void *p, size_t n_bytes)
{
	if(!p)
		return malloc(n_bytes);
	if(!n_bytes)
	{
		free(p);
		return NULL;
	}

	chunk_t *header = get_header(p);
	void *new_p = malloc(n_bytes);
	if(!new_p)
		return NULL;
	memcpy(new_p, p, header->size);
	free(p);
	return new_p;
}

size_t heap_usage(void)
{
	unsigned int usage = 0;
	chunk_t *head = (chunk_t *)&_heap_start;
	while(head->next != NULL)
	{
		if(!head->free)
			usage += head->size + sizeof(chunk_t);
		head = head->next;
	}
	return usage;
}
#endif
void *memcpy(void *i, const void *j, size_t size)
{
	const char *src = j;
	char *dst = i;
	while(size--)
	{
		*dst++ = *src++;
	}
	return i;
}

void *memset(void *i, int n, size_t size)
{
	char *p = i;
	while(size--)
	{
		*p++ = n;
	}
	return i;
}

size_t heap_usage(void)
{
	unsigned int usage = 0;
	chunk_t *head = (chunk_t *)&_heap_start;
	while(head->next != NULL)
	{
		if(!head->free)
			usage += head->size + sizeof(chunk_t);
		head = head->next;
	}
	return usage;
}

void __heap_status(int dump)
{
	size_t i;
	unsigned int cntr;
	unsigned int h_used, h_perc;
	h_used = heap_usage();
	h_perc = (h_used * 100)/(unsigned int)&_heap_size;
	printf("Heap Dump: %p - %p\n", &_heap_start, &_heap_end);
	printf("Heap Used: %u/%u - %u%%\n", h_used,
		(unsigned int)&_heap_size, h_perc);

	if(!dump)
		return;

	for(i = (size_t)&_heap_start; i < (size_t)&_heap_end; i+=32)
	{
		printf("[");
		for(cntr = 0; cntr < 32; cntr++)
		{
			printf("%02x", MMIO8(i+cntr));
			if(cntr && !((cntr + 1) % 4) && cntr != 31)
				printf(" ");
		}
		printf("]\n");
	}
}
int strncmp(const char *i, const char *j, size_t n)
{
	int ret = 0;
	if(!n)
		return 0;
	do
	{
		if(*i != *j++)
		{
			ret = (*(const unsigned char *)i - *(const unsigned char *)(j - 1));
			break;
		}
		if(*i++ == '\0')
		{
			ret = 0;
			break;
		}
	}
	while(--n);
	return ret;
}
static int __fputc(const int *dev, int en_stdout, const char c)
{
	int ret;
	putchar(c);
	if((c == '\n'))
		putchar('\r');
	return ret;
}

static int __fputs(const int *dev, int en_stdout, const char *i)
{
	int ret = 0;
	puts(i);
	return ret;
}
int _fputc(const FILE *dev, const char c)
{
	return __fputc(dev, 0, c);
}

int _fputs(const FILE *dev, const char *i)
{
	return __fputs(dev, 0, i);
}

static int unumprint(const int *dev, int en_stdout, unsigned long unum,
		unsigned int radix, char padc, int padn)
{
	char buf[20];
	int i = 0, ret = 0;
	do
	{
		unsigned int rem = unum % radix;
		if(rem < 0xa)
			buf[i] = '0' + (char)rem;
		else
			buf[i] = 'a' + (char)(rem - 0xa);
		i++;
		unum /= radix;
	}
	while (unum > 0U);
	if(padn > 0)
	{
		while(i < padn)
		{
			__fputc(dev, en_stdout, padc);
			ret++;
			padn--;
		}
	}
	while(--i >= 0)
	{
		__fputc(dev, en_stdout, buf[i]);
		ret++;
	}
	return ret;
}

#if USE_FLOAT == 1
static int fltprint(const FILE *dev, bool en_stdout, double flt,
		char padc, int padd, int padf)
{
	int ret = 0;
	long d = (long) flt;
	double frac = flt - (double) d;
	ret = unumprint(dev, en_stdout, d, 10, padc, padd);
	__fputc(dev, en_stdout, '.');
	ret ++;
	while(padf != 0)
	{
		frac *= 10.0;
		padf--;
	}
	d = (long) frac;
	ret += unumprint(dev, en_stdout, d,10, '0', 0);

	return ret;
}
#endif

int __vprintf(const int *dev, int en_stdout, const char *fmt, va_list args)
{
	int l_ret;
	long num;
	unsigned long unum;
	char *str;
	char ch;
	char padc = '\0';
	int padn;
	int ret = 0;
#if USE_FLOAT == 1
	double flt;
	int padf = 5;
#endif

	while(*fmt != '\0')
	{
		l_ret = 0;
		padn = 0;
		if(*fmt == '%')
		{
			fmt++;
loop:
			switch(*fmt)
			{
				case 'i':
				case 'd':
					num = get_num_va_args(args, l_ret);
					if (num < 0)
					{
						__fputc(dev, en_stdout, '-');
						unum = (unsigned long)-num;
						padn--;
					}
					else
						unum = (unsigned long)num;
					ret += unumprint(dev, en_stdout, unum, 10, padc, padn);
					break;
				case 'c':
					str = va_arg(args, char *);
					ret += __fputc(dev, en_stdout, (int)str);
					break;
				case 's':
					str = va_arg(args, char *);
					ret += __fputs(dev, en_stdout, str);
					break;
				case 'p':
					unum = (uintptr_t) va_arg(args, void *);
					ret += __fputs(dev, en_stdout, "0x");
					padn -= 2;
					ret += unumprint(dev, en_stdout, unum, 16, padc, padn);
					break;
				case 'x':
					unum = get_unum_va_args(args, l_ret);
					ret += unumprint(dev, en_stdout, unum, 16, padc, padn);
					break;
				case 'z':
					if (sizeof(size_t) == 8U)
						l_ret = 2;

					fmt++;
					goto loop;
				case 'l':
					l_ret++;
					fmt++;
					goto loop;
				case 'u':
					unum = get_unum_va_args(args, l_ret);
					ret += unumprint(dev, en_stdout, unum, 10, padc, padn);
					break;
#if USE_FLOAT == 1
				case 'f':
					flt = va_arg(args, double);
					if(flt < 0)
					{
						__fputc(dev, en_stdout, '-');
						flt *= -1.0;
						padn--;
					}
					ret += fltprint(dev, en_stdout, flt, padc, padn, padf);
					break;
#endif
				case '0':
					padc = '0';
					padn = 0;
					fmt++;
					while(1)
					{
						ch = *fmt;
						if((ch < '0') || (ch > '9'))
							goto loop;
						padn = (padn * 10) + (ch - '0');
						fmt++;
					}
#if USE_FLOAT == 1
				case '.':
					padf = 0;
					fmt++;
					while(true)
					{
						ch = *fmt;
						if((ch < '0') || (ch > '9'))
							goto loop;
						padf = (padf * 10) + (ch - '0');
						fmt++;
					}
#endif
				case '%':
					ret += __fputc(dev, en_stdout, *fmt);
					break;
				default:
					return -1;
			}
			fmt++;
			continue;
		}

		__fputc(dev, en_stdout, (char)*fmt);
		fmt++;
		ret++;
	}
	return ret;
}

int printf(const char *fmt, ...)
{
	int ret;
	va_list va;
	va_start(va, fmt);
	ret = __vprintf(0, 0, fmt, va);
	va_end(va);
	return ret;
}




#include <errno.h>
#ifndef _REENT_ONLY

int *
__errno ()
{
    int a=0;
  return &a;
}

#endif