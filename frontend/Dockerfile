# Stage 0, "build-stage", based on Node.js, to build and compile the frontend
FROM node:12.3.1-alpine as build-stage

WORKDIR /app

COPY package.json yarn.lock /app/

RUN yarn install --non-interactive

COPY ./ /app/

ARG FRONTEND_ENV=production

ENV VUE_APP_ENV=${FRONTEND_ENV}

# Comment out the next line to disable tests
RUN yarn test:unit

RUN yarn build



# Stage 1, based on Nginx, to have only the compiled app, ready for production with Nginx
FROM nginx:1.16.0-alpine

COPY --from=build-stage /app/dist/ /usr/share/nginx/html

COPY ./nginx.conf /etc/nginx/conf.d/default.conf
COPY ./nginx-backend-not-found.conf /etc/nginx/extra-conf.d/backend-not-found.conf
