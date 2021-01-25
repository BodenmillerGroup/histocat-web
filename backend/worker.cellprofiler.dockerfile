FROM plankter/histocat-worker-base:latest

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

ARG BACKEND_ENV=production

ENV BACKEND_ENV=${BACKEND_ENV}

COPY ./ /app
RUN chmod +x /app/bin/worker-start.sh

EXPOSE 5688

CMD ["/app/bin/worker-start.sh"]
 
