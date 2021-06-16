# Stage 0, "build-stage", based on Node.js, to build and compile the frontend
FROM node:14.17.1-alpine3.13 as build-stage

WORKDIR /app

COPY package.json yarn.lock /app/

RUN yarn install --frozen-lockfile --non-interactive

COPY ./ /app/

# Comment out the next line to disable tests
# RUN yarn test:unit

ARG VUE_APP_NAME
ENV VUE_APP_NAME=${VUE_APP_NAME}

RUN yarn build



# Stage 1, based on Nginx, to have only the compiled app, ready for production with Nginx
FROM nginx:1.20.0-alpine

LABEL maintainer="Anton Rau <anton.rau@uzh.ch>"

COPY --from=build-stage /app/dist/ /usr/share/nginx/html

COPY ./nginx.conf /etc/nginx/conf.d/default.conf
COPY ./nginx-backend-not-found.conf /etc/nginx/extra-conf.d/backend-not-found.conf
