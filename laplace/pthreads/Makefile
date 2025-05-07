CC = gcc
CFLAGS = -pthread

all: clean iterative pthread

iterative:
	@cd sequential && \
	$(CC) seq.c -o s_laplace && \
	$(CC) seq_lin.c -o s_laplace_lin && \
	$(CC) seq_ptr.c -o s_laplace_ptr && \
	$(CC) seq_linptr.c -o s_laplace_linptr
	@mv sequential/s_* binaries

pthread:

	@cd demand && \
	$(CC) $(CFLAGS) pthreads.c -o d_laplace && \
	$(CC) $(CFLAGS) pthread_lin.c -o d_laplace_lin && \
	$(CC) $(CFLAGS) pthread_ptr.c -o d_laplace_ptr && \
	$(CC) $(CFLAGS) pthread_linptr.c -o d_laplace_linptr
	@mv demand/d_* binaries

	@cd regular && \
	$(CC) $(CFLAGS) pthreads.c -o r_laplace && \
	$(CC) $(CFLAGS) pthread_lin.c -o r_laplace_lin && \
	$(CC) $(CFLAGS) pthread_ptr.c -o r_laplace_ptr && \
	$(CC) $(CFLAGS) pthread_linptr.c -o r_laplace_linptr
	@mv regular/r_* binaries

	@cd round_robin && \
	$(CC) $(CFLAGS) pthreads.c -o rr_laplace && \
	$(CC) $(CFLAGS) pthread_lin.c -o rr_laplace_lin && \
	$(CC) $(CFLAGS) pthread_ptr.c -o rr_laplace_ptr && \
	$(CC) $(CFLAGS) pthread_linptr.c -o rr_laplace_linptr
	@mv round_robin/rr_* binaries

clean:
	@find . -type f -executable -exec rm -f {} +
