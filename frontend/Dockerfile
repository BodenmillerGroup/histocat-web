# Stage 0, "build-stage", based on Node.js, to build and compile the frontend
FROM node:14.15.1-alpine3.12 as build-stage

# Comment out the next line to disable git
#RUN apk update && apk upgrade && \
#    apk add --no-cache git

WORKDIR /app

COPY package.json yarn.lock /app/

RUN yarn install --non-interactive

COPY ./ /app/

ARG FRONTEND_ENV
ENV VUE_APP_ENV=${FRONTEND_ENV}

# Comment out the next line to disable tests
RUN yarn test:unit

RUN yarn build



# Stage 1, based on Nginx, to have only the compiled app, ready for production with Nginx
FROM nginx:1.18.0-alpine

LABEL maintainer="Anton Rau <anton.rau@gmail.com>"

COPY --from=build-stage /app/dist/ /usr/share/nginx/html

COPY ./nginx.conf /etc/nginx/conf.d/default.conf
COPY ./nginx-backend-not-found.conf /etc/nginx/extra-conf.d/backend-not-found.conf
