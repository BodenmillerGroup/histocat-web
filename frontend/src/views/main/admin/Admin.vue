<template>
  <router-view></router-view>
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import { store } from "@/store";
import { Component, Vue } from "vue-property-decorator";

const mainContext = mainModule.context(store);

const routeGuardAdmin = async (to, from, next) => {
  if (!mainContext.getters.hasAdminAccess) {
    next("/main");
  } else {
    next();
  }
};

@Component
export default class Admin extends Vue {
  beforeRouteEnter(to, from, next) {
    routeGuardAdmin(to, from, next);
  }

  beforeRouteUpdate(to, from, next) {
    routeGuardAdmin(to, from, next);
  }
}
</script>
