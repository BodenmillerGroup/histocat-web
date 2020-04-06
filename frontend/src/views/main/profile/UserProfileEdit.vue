<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Edit User Profile</div>
      </v-card-title>
      <v-card-text>
        <template>
          <v-form v-model="valid" ref="form" lazy-validation>
            <v-text-field label="Name" v-model="name"></v-text-field>
            <v-text-field label="E-mail" type="email" v-model="email" :rules="emailRules"></v-text-field>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { mainModule } from "@/modules/main";
import { required, email } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";
import { IUserProfileUpdate } from "@/modules/user/models";

@Component
export default class UserProfileEdit extends Vue {
  readonly mainContext = mainModule.context(this.$store);

  readonly emailRules = [required, email];

  valid = true;
  name: string = "";
  email: string = "";

  created() {
    const userProfile = this.userProfile;
    if (userProfile) {
      this.name = userProfile.name;
      this.email = userProfile.email;
    }
  }

  get userProfile() {
    return this.mainContext.getters.userProfile;
  }

  reset() {
    const userProfile = this.userProfile;
    if (userProfile) {
      this.name = userProfile.name;
      this.email = userProfile.email;
    }
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const updatedProfile: IUserProfileUpdate = {};
      if (this.name) {
        updatedProfile.name = this.name;
      }
      if (this.email) {
        updatedProfile.email = this.email;
      }
      await this.mainContext.actions.updateUserProfile(updatedProfile);
      this.$router.push("/main/profile");
    }
  }
}
</script>
