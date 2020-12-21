/**
 * Welcome to your Workbox-powered service worker!
 *
 * You'll need to register this file in your web app and you should
 * disable HTTP caching for this file too.
 * See https://goo.gl/nhQhGp
 */

//Change this value every time before you build
const LATEST_VERSION = "2020.12.21";

workbox.core.setCacheNameDetails({ prefix: "histocat" });

self.addEventListener("activate", (event) => {
  console.log(`%c ${LATEST_VERSION} `, "background: #ddd; color: #0000ff");
  if (caches) {
    caches.keys().then((arr) => {
      arr.forEach((key) => {
        if (key.indexOf("histocat-precache") < -1) {
          caches.delete(key).then(() => console.log(`%c Cleared ${key}`, "background: #333; color: #ff0000"));
        } else {
          caches.open(key).then((cache) => {
            cache.match("version").then((res) => {
              if (!res) {
                cache.put("version", new Response(LATEST_VERSION, { status: 200, statusText: LATEST_VERSION }));
              } else if (res.statusText !== LATEST_VERSION) {
                caches
                  .delete(key)
                  .then(() => console.log(`%c Cleared Cache ${LATEST_VERSION}`, "background: #333; color: #ff0000"));
              } else
                console.log(`%c You have the latest version ${LATEST_VERSION}`, "background: #333; color: #00ff00");
            });
          });
        }
      });
    });
  }
});

workbox.skipWaiting();
workbox.clientsClaim();

/**
 * The workboxSW.precacheAndRoute() method efficiently caches and responds to
 * requests for URLs in the manifest.
 * See https://goo.gl/S9QRab
 */
self.__precacheManifest = [].concat(self.__precacheManifest || []);
// workbox.precaching.suppressWarnings();
workbox.precaching.precacheAndRoute(self.__precacheManifest, {});
