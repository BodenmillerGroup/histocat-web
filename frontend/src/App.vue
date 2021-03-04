<template>
  <v-app id="app">
    <v-main v-if="loggedIn === null">
      <v-container fill-height>
        <v-row align="center" justify="center">
          <v-col>
            <div class="text-center">
              <div class="text-h5 my-12">Loading...</div>
              <v-progress-circular size="100" indeterminate color="primary"></v-progress-circular>
            </div>
          </v-col>
        </v-row>
      </v-container>
    </v-main>
    <router-view v-else />
    <NotificationsManager />
  </v-app>
</template>

<script lang="ts">
import NotificationsManager from "@/components/NotificationsManager.vue";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";
import { responsiveModule } from "@/modules/responsive";

@Component({
  components: {
    NotificationsManager,
  },
})
export default class App extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly responsiveContext = responsiveModule.context(this.$store);

  get loggedIn() {
    return this.mainContext.getters.isLoggedIn;
  }

  async created() {
    await this.mainContext.actions.checkLoggedIn();
  }

  mounted() {
    /* listen for resize events */
    window.addEventListener("resize", () => {
      this.responsiveContext.mutations.setResponsive({
        height: window.innerHeight,
        width: window.innerWidth,
      });
    });
    this.responsiveContext.mutations.setResponsive({
      height: window.innerHeight,
      width: window.innerWidth,
    });
  }
}
</script>
