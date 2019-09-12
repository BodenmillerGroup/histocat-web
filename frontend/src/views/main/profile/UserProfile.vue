<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">User Profile</div>
      </v-card-title>
      <v-card-text>
        <div class="my-6">
          <div class="subtitle-1 primary--text text--lighten-3">Full Name</div>
          <div class="title primary--text text--darken-2" v-if="userProfile && userProfile.full_name">
            {{ userProfile.full_name }}
          </div>
          <div class="title primary--text text--darken-2" v-else>-----</div>
        </div>
        <div class="my-4">
          <div class="subtitle-1 primary--text text--lighten-3">Email</div>
          <div class="title primary--text text--darken-2" v-if="userProfile && userProfile.email">
            {{ userProfile.email }}
          </div>
          <div class="title primary--text text--darken-2" v-else>-----</div>
        </div>
      </v-card-text>
      <v-card-actions>
        <v-btn to="/main/profile/edit">Edit</v-btn>
        <v-btn to="/main/profile/password">Change password</v-btn>
        <v-spacer></v-spacer>
        <v-btn color="secondary" @click="resetSettings">Reset Settings</v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import { settingsModule } from "@/modules/settings";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class UserProfile extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly settingsModule = settingsModule.context(this.$store);

  get userProfile() {
    return this.mainContext.getters.userProfile;
  }

  resetSettings() {
    if (self.confirm("Reset all settings?")) {
      this.settingsModule.actions.resetSettings();
    }
  }
}
</script>
