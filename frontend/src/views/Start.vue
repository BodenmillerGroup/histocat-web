<template>
  <router-view />
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import { store } from "@/store";
import { Component, Vue } from "vue-property-decorator";

const mainContext = mainModule.context(store);

const startRouteGuard = async (to, from, next) => {
  await mainContext.actions.checkLoggedIn();
  if (mainContext.getters.isLoggedIn) {
    if (to.path === "/login" || to.path === "/") {
      next("/main");
    } else {
      next();
    }
  } else if (mainContext.getters.isLoggedIn === false) {
    if (to.path === "/" || (to.path as string).startsWith("/main")) {
      next("/login");
    } else {
      next();
    }
  }
};

@Component
export default class Start extends Vue {
  beforeRouteEnter(to, from, next) {
    startRouteGuard(to, from, next);
  }

  beforeRouteUpdate(to, from, next) {
    startRouteGuard(to, from, next);
  }
}
</script>
